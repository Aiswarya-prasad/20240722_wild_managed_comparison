"""
name: backmapping-binning
description:
"""

"""
Pre-backmapping:
Use simka to compare k-mer similarity between samples in order to choose 50 samples to map to each assembly.
"""

# subset with simka
# compare k-mer similarity between samples in order to choose 50 samples to map to each assembly.

rule make_simka_input_list:
    input:
        reads1 = expand("results/02_HostDepleting/{sample}/{sample}_R1_filt.fq.gz", sample=SAMPLES),
        reads2 = expand("results/02_HostDepleting/{sample}/{sample}_R2_filt.fq.gz", sample=SAMPLES),
    output:
        simka_list = "results/04_binning_backmapping/simka_input_list.txt"
    log: "results/04_binning_backmapping/simka_input_list.log"
    # conda: "../config/envs/scripts-env.yaml"
    run:
        '''
        Create a list of reads files for simka
        ID1: filename_pair1.fasta ; filename_pair2.fasta
        where ID is the sample name
        '''
        with open(output.simka_list, "w") as f:
            for r1, r2 in zip(input.reads1, input.reads2):
                sample = os.path.basename(r1).split("_R")[0]
                f.write(f"{sample}: {os.path.join(PROJECT_FULL_PATH, r1)}; {os.path.join(PROJECT_FULL_PATH, r2)}\n")

rule simka:
    input:
        simka_list = "results/04_binning_backmapping/simka_input_list.txt"
    output:
        simka_out = directory("results/04_binning_backmapping/simka_out")
    params:
        jobname="simka",
        account="pengel_spirit",
        runtime_s=convertToSec("0-23:10:00"),
        tmpdir = directory("/scratch/aprasad/20240722_wild_managed_comparison_tmp/simka_tmp"),
        maxreads = 0,
        abun = 2,
        maxcount = 200,
        maxmerge = 200,
        mem = 19000
    resources:
        mem_mb = convertToMb("20G")
    threads: 8
    log: "results/04_binning_backmapping/simka.log"
    conda: "../config/envs/mags-env.yaml" # actually used slightly different "../config/envs/mags-env.yaml"
    shell:
        """
        simka -in {input.simka_list} -out {output.simka_out} \
            -out-tmp {params.tmpdir} -keep-tmp \
            -max-reads {params.maxreads} -abundance-min {params.abun} \
            -max-count {params.maxcount} -max-merge {params.maxmerge} \
            -max-memory {params.mem} -nb-cores {threads} &> {log}
        """

# This rule and script are from Meline's pipeline: https://github.com/momsane/Metabee

rule parse_simka:
    input:
        simka_csv = "results/04_binning_backmapping/simka_out/mat_abundance_{dist}.csv.gz"
    output:
        heatmap = "results/04_binning_backmapping/simka_parsed/heatmap_{dist}.png",
        table = "results/04_binning_backmapping/simka_parsed/top_50_similar_{dist}.tsv"
    params:
        jobname="parse_simka",
        account="pengel_spirit",
        runtime_s=convertToSec("0-23:10:00"),
    resources:
        mem_mb = convertToMb("10G")
    log: "results/04_binning_backmapping/simka_parsed/parse_simka_{dist}.log"
    conda: "../config/envs/rmd-env.yaml"
    shell:
        """
        Rscript --vanilla scripts/simka_heatmap_similar.R {input.simka_csv} {output.heatmap} {output.table} {wildcards.dist}
        """

"""
Backmapping:
Map reads from selected samples to each assembly.
"""

rule index_bwa:
    input:
        assembly = "results/03_Assembly/{sample}/{sample}_scaffolds_1k_filtered.fasta"
    output:
        index = multiext("results/03_Assembly/{sample}/{sample}_scaffolds_1k_filtered.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
    params:
        jobname="index_bowtie_{sample}",
        account="pengel_spirit",
        runtime_s=convertToSec("0-02:00:00"),
    resources:
        mem_mb = convertToMb("10G")
    threads: 4
    log: "results/03_Assembly/{sample}/{sample}_scaffolds_1k_filtered_bowtie_index.log"
    benchmark: "results/03_Assembly/{sample}/{sample}_scaffolds_1k_filtered_bowtie_index.benchmark"
    conda: "../config/envs/mapping-bowtie-env.yaml"
    shell:
        """
        bwa index {input.assembly} &> {log}
        """

# backmapping on selected subset

rule backmapping:
    priority: 10
    input:
        assembly = "results/03_Assembly/{assembly}/{assembly}_scaffolds_1k_filtered.fasta",
        index = multiext("results/03_Assembly/{sample}/{sample}_scaffolds_1k_filtered.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
        reads1 = "results/02_HostDepleting/{sample}/{sample}_R1_filt.fq.gz",
        reads2 = "results/02_HostDepleting/{sample}/{sample}_R2_filt.fq.gz"
    output:
        bam = temp(SCRATCH_PATH+"/backmapping_tmp/{assembly}/{assembly}_scaffolds_vs_{sample}.bam")
    params:
        jobname="backmapping_{sample}",
        account="pengel_spirit",
        runtime_s=convertToSec("0-15:10:00"),
        tmpdir_prefix = SCRATCH_PATH+"/backmapping_tmp/tmp",
    resources:
        mem_mb = convertToMb("50G")
    threads: 16
    log: "results/04_binning_backmapping/backmapping/{assembly}/{assembly}_scaffolds_vs_{sample}.log"
    benchmark: "results/04_binning_backmapping/backmapping/{assembly}/{assembly}_scaffolds_vs_{sample}.benchmark"
    conda: "../config/envs/mapping-bowtie-env.yaml"
    shell:
        """
        mkdir -p {params.tmpdir_prefix}
        bwa mem -a -t {threads} {input.assembly} \
            {input.reads1} {input.reads2} | \
            samtools view -F 4 -bh - | \
            samtools sort -T {params.tmpdir_prefix} -O bam -@ {threads} > {output.bam}
        """

rule make_depth_files:
    priority: 20
    input:
        bam = SCRATCH_PATH+"/backmapping_tmp/{assembly}/{assembly}_scaffolds_vs_{sample}.bam"
    output:
        depth_file = SCRATCH_PATH+"/backmapping_tmp/{assembly}/{assembly}_scaffolds_vs_{sample}.depth"
    params:
        jobname="depth_{sample}",
        account="pengel_spirit",
        runtime_s=convertToSec("0-20:10:00"),
        tmpdir_prefix = SCRATCH_PATH+"/backmapping_tmp/tmp"
    resources:
        mem_mb = convertToMb("50G")
    threads: 16
    log: "results/04_binning_backmapping/backmapping/{assembly}/{assembly}_scaffolds_vs_{sample}_depth_summary.log"
    benchmark: "results/04_binning_backmapping/backmapping/{assembly}/{assembly}_scaffolds_vs_{sample}_depth_summary.benchmark"
    conda: "../config/envs/mags-env.yaml" # actually used slightly different "../config/envs/mags-env.yaml"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output.depth_file} {input.bam} &> {log}
        """

rule merge_depth_files:
    priority: 30
    input:
        combinations = "results/04_binning_backmapping/simka_parsed/top_50_similar_jaccard.tsv",
        depth_files = lambda wildcards: top50_samples_for_backmapping(comb_file = "results/04_binning_backmapping/simka_parsed/top_50_similar_jaccard.tsv",
                                                                      depth_files_prefix = SCRATCH_PATH+"/backmapping_tmp/{assembly}/{assembly}_scaffolds_vs_",
                                                                      assembly_name = wildcards.assembly),
    output:
        merged_depth_file = "results/04_binning_backmapping/backmapping/{assembly}/{assembly}_scaffolds_vs_top50_samples.depth"
    params:
        jobname="merge_depth_{assembly}",
        account="pengel_spirit",
        runtime_s=convertToSec("0-15:10:00"),
        tmpdir_prefix = SCRATCH_PATH+"/backmapping_tmp/tmp"
    resources:
        mem_mb = convertToMb("50G")
    threads: 16
    log: "results/04_binning_backmapping/backmapping/{assembly}/{assembly}_scaffolds_vs_all_samples_depth_summary.log"
    benchmark: "results/04_binning_backmapping/backmapping/{assembly}/{assembly}_scaffolds_vs_all_samples_depth_summary.benchmark"
    conda: "../config/envs/mags-env.yaml" # actually used slightly different "../config/envs/mags-env.yaml"
    shell:
        """
        ./scripts/merge_depths.pl {input.depth_files} > {output.merged_depth_file}
        """

# binning on all samples with metabat
rule bin_metabat:
    input:
        merged_depth_file = "results/04_binning_backmapping/backmapping/{sample}/{sample}_scaffolds_vs_top50_samples.depth",
        assembly = "results/03_Assembly/{sample}/{sample}_scaffolds_1k_filtered.fasta",
    output:
        out_dir = directory("results/04_binning_backmapping/binning/{sample}_MAGs")
    params:
        min_len = 2000,
        max_edges = 500,
        min_cls = 200000,
        min_cv = 1,
        jobname="binning_{sample}",
        account="pengel_spirit",
        runtime_s=convertToSec("0-10:10:00"),
    resources:
        mem_mb = convertToMb("50G")
    threads: 16
    log: "results/04_binning_backmapping/binning/{sample}_mag_binning.log"
    benchmark: "results/04_binning_backmapping/binning/{sample}_mag_binning.benchmark"
    conda: "../config/envs/mags-env.yaml" # actually used slightly different "../config/envs/mags-env.yaml"
    shell:
        """
        metabat2 -i {input.assembly} -a {input.merged_depth_file} \
        -o {output.out_dir}/{wildcards.sample} --minContig {params.min_len} \
        --maxEdges {params.max_edges} -x {params.min_cv} --minClsSize {params.min_cls} --saveCls -v
        """