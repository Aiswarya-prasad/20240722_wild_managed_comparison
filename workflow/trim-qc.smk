#!/usr/bin/env python

"""
name: trim-qc
description: 
"""

rule raw_qc:
    input:
        reads = lambda wildcards: get_input_file(sample = wildcards.sample, read = wildcards.read, lane = wildcards.lane),
    output:
        html = "results/00_RawQC/{sample}_{lane}_{read}_fastqc.html",
        zip = "results/00_RawQC/{sample}_{lane}_{read}_fastqc.zip"
    params:
        outdir = lambda wildcards: f"results/00_RawQC/{wildcards.sample}",
        outname = lambda wildcards: f"{wildcards.sample}_{wildcards.lane}_{wildcards.read}",
        mailto = "aiswarya.prasad@unil.ch",
        mailtype = "BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname = lambda wildcards: f"{wildcards.sample}_{wildcards.lane}_{wildcards.read}_qc",
        account = "pengel_spirit",
        runtime_s = convertToSec("0-2:10:00"),
    resources:
        mem_mb = 8000
    log: "results/00_RawQC/{sample}/{sample}_{lane}_{read}_qc.log"
    benchmark: "results/00_RawQC/{sample}/{sample}_{lane}_{read}_qc.benchmark"
    threads: 2
    conda: "../config/envs/trim-qc-env.yaml"
    shell:
        """
        zcat {input.reads} | fastqc -t {threads} -o {params.outdir} stdin:{params.outname}
        """

"""
information from novogene:
If we agree to deliver the clean data before the project starts, we will filter the data strictly according to the standard to obtain high quality clean data which can be used for further research and paper writing. We will discard the paired reads in the following situation: when either one read contains adapter contamination; when either one read contains more than 10 percent uncertain nucleotides; when either one read contains more than 50 percent low quality nucleotides (base quality less than 5). The data analysis results based on the clean data that is filtered by this standard can be approved by high level magazines (Yan L.Y. et al . 2013). If you want to get more information, please refer to the official website of Novogene (www.novogene.com).

1.4 Results of Raw Data Filtering
The sequenced reads (raw reads) often contain low quality reads and adapters, which will affect the analysis quality. So it's necessary to filter the raw reads and get the clean reads. The filtering process is as follows:

(1) Remove reads containing adapters.

(2) Remove reads containing N > 10% (N represents the base cannot be determined).

(3) Remove reads containing low quality (Qscore<= 5) base which is over 50% of the total base.

    Sequences of adapter

    5' Adapter:

    5'-AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT-3'

    3' Adapter:

    5'-GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG-3'
"""

rule trim:
    input:
        reads1 = lambda wildcards: get_input_file(sample = wildcards.sample, read = "1", lane = wildcards.lane),
        reads2 = lambda wildcards: get_input_file(sample = wildcards.sample, read = "2", lane = wildcards.lane),
        adapters = "data/AdaptersPE-mod.fa"
    output:
        reads1 = "results/01_TrimmingFiltering/{sample}_{lane}_1.fq.gz",
        reads2 = "results/01_TrimmingFiltering/{sample}_{lane}_2.fq.gz",
        reads1_unpaired = temp("results/01_TrimmingFiltering/{sample}_{lane}_1_unpaired.fq.gz"),
        reads2_unpaired = temp("results/01_TrimmingFiltering/{sample}_{lane}_2_unpaired.fq.gz")
    params:
        jobname="{sample}_{lane}_trim",
        account="pengel_spirit",
        runtime_s=convertToSec("0-4:10:00"),
    resources:
        mem_mb = 8000
    threads: 4
    log: "results/00_trimmedreads/{sample}_{lane}_trim.log"
    benchmark: "results/00_trimmedreads/{sample}_{lane}_trim.benchmark"
    conda: "../config/envs/trim-qc-env.yaml"
    shell:
        """
        trimmomatic PE -threads {threads} {input.reads1} {input.reads2} \
            {output.reads1} {output.reads1_unpaired} \
            {output.reads2} {output.reads2_unpaired} \
            ILLUMINACLIP:{input.adapters}:2:30:10 LEADING:28 \
            TRAILING:28  MINLEN:60 &> {log}
        """

rule trim_qc:
    input:
        reads = "results/01_TrimmingFiltering/{sample}_{lane}_{read}.fq.gz",
    output:
        html="results/01_TrimmedQC/{sample}_{lane}_{read}_fastqc.html",
        zip="results/01_TrimmedQC/{sample}_{lane}_{read}_fastqc.zip"
    threads: 2
    conda: "../config/envs/trim-qc-env.yaml"
    params:
        outdir="results/01_TrimmedQC/",
        jobname="{sample}_{lane}_{read}_trim_qc",
        account="pengel_spirit",
        runtime_s=convertToSec("0-4:10:00"),
    log: "results/01_TrimmingFiltering/fastqc/{sample}/{sample}_{lane}_{read}_trim_qc.log"
    benchmark: "results/01_TrimmingFiltering/fastqc/{sample}/{sample}_{lane}_{read}_trim_qc.benchmark"
    resources:
        mem_mb = 8000
    shell:
        """
        fastqc -t {threads} {input.reads} -o {params.outdir} &> {log}
        """

rule get_fastqc_text_outputs_raw:
    input:
        fastqc_zip = "results/00_RawQC/{sample}_{lane}_{read}_fastqc.zip"
    output:
        text = "results/00_RawQC/{sample}_{lane}_{read}_fastqc.txt"
    params:
        jobname="{sample}_{lane}_{read}_get_fastqc_text",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
    resources:
        mem_mb = 8000
    log: "results/00_RawQC/{sample}/{sample}_{lane}_{read}_get_fastqc_text.log"
    benchmark: "results/00_RawQC/{sample}/{sample}_{lane}_{read}_get_fastqc_text.benchmark"
    shell:
        """
        unzip -p {input.fastqc_zip} */fastqc_data.txt > {output.text}
        """

rule get_fastqc_text_outputs_trimmed:
    input:
        fastqc_zip = "results/01_TrimmedQC/{sample}_{lane}_{read}_fastqc.zip"
    output:
        text = "results/01_TrimmedQC/{sample}_{lane}_{read}_fastqc.txt"
    params:
        jobname="{sample}_{lane}_{read}_get_fastqc_text",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
    resources:
        mem_mb = 8000
    log: "results/01_TrimmedQC/{sample}/{sample}_{lane}_{read}_get_fastqc_text.log"
    benchmark: "results/01_TrimmedQC/{sample}/{sample}_{lane}_{read}_get_fastqc_text.benchmark"
    shell:
        """
        unzip -p {input.fastqc_zip} */fastqc_data.txt > {output.text}
        """

rule fetch_host_genome:
    output:
        host_genome = "data/host_genome/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
    params:
        host_genome_link = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz",
        jobname="fetch_host_genome",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
    resources:
        mem_mb = 8000
    log: "data/host_genome/fetch_host_genome.log"
    benchmark: "data/host_genome/fetch_host_genome.benchmark"
    shell:
        """
        wget {params.host_genome_link} -O {output.host_genome}.gz &> {log}
        echo "Host genome downloaded"
        gunzip {output.host_genome}.gz &> {log}
        """


rule bowtie_index:
    input:
        host_database = "data/host_genome/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
    output:
        bowtie_index = multiext("data/host_genome/GCF_003254395.2_Amel_HAv3.1_genomic.fna", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    params:
        jobname="bowtie_index",
        account="pengel_spirit",
        runtime_s=convertToSec("0-07:10:00"),
    resources:
        mem_mb = convertToMb("50G")
    threads: 4
    log: "results/02_HostFiltering/bowtie_index.log"
    benchmark: "results/02_HostFiltering/bowtie_index.benchmark"
    conda: "../config/envs/mapping-bowtie-env.yaml"
    shell:
        """
        bowtie2-build {input.host_database} {input.host_database} &> {log}
        """

rule map_to_host_bowtie2:
    input:
        reads1 = lambda wildcards: get_trimmed_files(sample = wildcards.sample, read = "1"),
        reads2 = lambda wildcards: get_trimmed_files(sample = wildcards.sample, read = "2"),
        bowtie_index = multiext("data/host_genome/GCF_003254395.2_Amel_HAv3.1_genomic.fna", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
        host_database = "data/host_genome/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
    output:
        bam = temp("results/02_HostFiltering/{sample}/{sample}_bowtie.bam"),
        flagstat = "results/02_HostFiltering/{sample}/{sample}_bowtie_flagstat.tsv",
    params:
        jobname="{sample}_map_to_host_bowtie",
        account="pengel_spirit",
        runtime_s=convertToSec("0-17:10:00"),
        input_files_string1 = lambda wildcards: ','.join(get_trimmed_files(sample = wildcards.sample, read = "1")),
        input_files_string2 = lambda wildcards: ','.join(get_trimmed_files(sample = wildcards.sample, read = "2")),
    resources:
        mem_mb = convertToMb("50G")
    threads: 4
    log: "results/02_HostFiltering/{sample}/{sample}_map_to_rep_MAGs_bowtie.log"
    benchmark: "results/02_HostFiltering/{sample}/{sample}_map_to_rep_MAGs_bowtie.benchmark"
    conda: "../config/envs/mapping-bowtie-env.yaml"
    shell:
        """
        bowtie2 -x {input.host_database} -1 {params.input_files_string1} \
            -2 {params.input_files_string1} \
            | samtools view -bh - | samtools sort - > {output.bam}
        samtools flagstat {output.bam} > {output.flagstat}
        """

rule get_non_host_reads:
    input:
        bam = "results/02_HostFiltering/{sample}/{sample}_bowtie.bam",
    output:
        reads1 = "results/02_HostFiltering/{sample}/{sample}_R1_hostfilt.fq.gz",
        reads2 = "results/02_HostFiltering/{sample}/{sample}_R2_hostfilt.fq.gz",
        bam_unmapped = temp("results/02_HostFiltering/{sample}/{sample}_non_host_unmapped.bam"),
    params:
        jobname="{sample}_get_non_host_reads",
        account="pengel_spirit",
        runtime_s=convertToSec("0-07:10:00"),
    resources:
        mem_mb = convertToMb("50G")
    threads: 4
    log: "results/02_HostFiltering/{sample}/{sample}_get_non_host_reads.log"
    benchmark: "results/02_HostFiltering/{sample}/{sample}_get_non_host_reads.benchmark"
    conda: "../config/envs/mapping-bowtie-env.yaml"
    shell:
        """
        samtools view -b -f 4 {input.bam} | samtools sort -n - > {output.bam_unmapped}
        outreads1={output.reads1}
        outreads2={output.reads2}
        bedtools bamtofastq -i {input.bam} -fq ${{outreads1/.gz/}} -fq2 ${{outreads2/.gz/}} &> {log}
        gzip ${{outreads1/.gz/}}
        gzip ${{outreads2/.gz/}}
        """

# rule concatenate_reads:
#     input:
#         reads = lambda wildcards: [f"results/00_trimmedreads/{x.split('.fq.gz')[0]}_trim.fq.gz" for x in get_list_of_values(get_renamed_input_files({key: raw_paths_dict[key] for key in raw_paths_dict.keys() if key == wildcards.sample})) if f"_{wildcards.read}" in x]
#     output:
#         concat_reads = "results/01_trimmedconcatreads/{sample}_{read}.fq.gz"
#     threads: 2
#     conda: "../config/envs/trim-qc-env.yaml"
#     params:
#      
#         jobname="{sample}_{read}_concat",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-2:10:00"),
#     log: "results/01_trimmedconcatreads/{sample}_{read}_concat.log"
#     benchmark: "results/01_trimmedconcatreads/{sample}_{read}_concat.benchmark"
#     resources:
#         mem_mb = 8000
#     shell:
#         """
#         cat {input.reads} > {output.concat_reads}
#         """


# rest of the pipeline starts from clean reads so commenting out this rule for now
# rule re_pair_reads:
#     input:
#         reads1 = "results/01_trimmedconcatreads/{sample}_R1.fq.gz",
#         reads2 = "results/01_trimmedconcatreads/{sample}_R2.fq.gz"
#     output:
#         reads1 = "results/01_cleanreads/{sample}_R1_repaired.fq.gz",
#         reads2 = "results/01_cleanreads/{sample}_R2_repaired.fq.gz",
#         singletons = "results/01_cleanreads/{sample}_singletons.fq.gz"
#     params:
#         java_mem="170", # 70 worked for most but 11/150 samples needed more - should not be more than 85% total memory requested
#      
#         jobname="{sample}_re-paired",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-23:00:00"),
#     resources:
#         mem_mb = convertToMb("400G")
#     threads: 2 # 4 worked for most but 11/150 samples needed more
#     log: "results/01_cleanreads/{sample}_repaired.log"
#     benchmark: "results/01_cleanreads/{sample}_repaired.benchmark"
#     conda: "../config/envs/mapping-bowtie-env.yaml"
#     shell:
#         """
#         repair.sh -Xmx{params.java_mem}g threads={threads} in1={input.reads1} in2={input.reads2} out1={output.reads1} out2={output.reads2} outs={output.singletons} repair &> {log}
#         """