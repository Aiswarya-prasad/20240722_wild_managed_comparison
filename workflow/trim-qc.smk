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