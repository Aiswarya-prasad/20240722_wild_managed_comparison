"""
name: backmapping-binning
description:
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
    conda: "../config/envs/mags-env.yaml"
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
        heatmap = "results/04_binning_backmapping/simka_parsed/heatmap_{dist}.pdf"
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

# This rule is from Meline's pipeline: https://github.com/momsane/Metabee

rule make_combinations:
    input:
        table = "results/04_binning_backmapping/simka_parsed/top_50_similar_{dist}.tsv"
    output:
        result = "results/04_binning_backmapping/simka_parsed/combinations_{dist}.txt"
    threads: 1
    resources:
        account = "pengel_general_data",
        runtime = "20m",
        mem_mb = 500
    log: "results/04_binning_backmapping/simka_parsed/combinations_{dist}.log"
    shell:
        """
        cat {input.table} | awk ' NR > 1 {{ print "A_"$1"_R_"$2 }}' > {output.results};
        """

# backmapping on selected subset

# binning on all samples with metabat