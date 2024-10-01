#!/usr/bin/env python

"""
name: assembly_qc
description: 
"""

rule run_bbnorm:
    input:
        reads1 = "results/02_HostDepleting/{sample}/{sample}_R1_filt.fq.gz",
        reads2 = "results/02_HostDepleting/{sample}/{sample}_R2_filt.fq.gz",
    output:
        reads1_norm = "results/03_Assembly/{sample}/{sample}_R1_filt_norm.fq.gz",
        reads2_norm = "results/03_Assembly/{sample}/{sample}_R2_filt_norm.fq.gz",
        hist = "results/03_Assembly/{sample}/{sample}_bbnorm_hist.txt",
        peaks = "results/03_Assembly/{sample}/{sample}_bbnorm_peaks.txt",
    params:
        jobname="{sample}_bbnorm",
        account="pengel_spirit",
        runtime_s=convertToSec("0-23:10:00"), # increased for sample DNA27_2 which did not finish in 7 hours
        # runtime_s=convertToSec("0-07:10:00"),
        java_mem=50,
    resources:
        mem_mb = convertToMb("250G")
    threads: 8 # increased for sample DNA27_2 which did not finish in 23 hours with 4 threads
    log: "results/03_Assembly/{sample}/{sample}_bbnorm.log"
    benchmark: "results/03_Assembly/{sample}/{sample}_bbnorm.benchmark"
    conda: "../config/envs/mapping-bowtie-env.yaml"
    shell:
        """
        bbnorm.sh -Xmx{params.java_mem}G threads={threads} hist={output.hist} peaks={output.peaks} \
            target=40 mindepth=0 \
            in1={input.reads1} in2={input.reads2} out1={output.reads1_norm} out2={output.reads2_norm} &> {log}
        """

rule run_spades:
    input:
        reads1 = "results/03_Assembly/{sample}/{sample}_R1_filt_norm.fq.gz",
        reads2 = "results/03_Assembly/{sample}/{sample}_R2_filt_norm.fq.gz",
    output:
        contigs = "results/03_Assembly/{sample}/{sample}_contigs.fasta",
        scaffolds = "results/03_Assembly/{sample}/{sample}_scaffolds.fasta",
        graph = "results/03_Assembly/{sample}/{sample}_assembly_graph.fastg",
        spades_log = "results/03_Assembly/{sample}/{sample}_spades.log",
    params:
        jobname="{sample}_spades",
        account="pengel_spirit",
        runtime_s=convertToSec("0-23:10:00"),
        outdir = "results/03_Assembly/{sample}_assembly_dir/",
    resources:
        spades_mem=lambda wildcards, attempt: f'{250 * attempt}',
        mem_mb = lambda wildcards, attempt: convertToMb(f'{250 * attempt}G'),
    retries: 2
    threads: 4
    log: "results/03_Assembly/{sample}/{sample}_run_spades.log"
    benchmark: "results/03_Assembly/{sample}/{sample}_run_spades.benchmark"
    conda: "../config/envs/spades-env-updated.yaml"
    shell:
        """
        spades.py -t {threads} -m {resources.spades_mem} --only-assembler \
            --pe1-1 {input.reads1} --pe1-2 {input.reads2} \
            -o {params.outdir}
        mv results/03_Assembly/{wildcards.sample}_assembly_dir/contigs.fasta {output.contigs}
        mv results/03_Assembly/{wildcards.sample}_assembly_dir/scaffolds.fasta {output.scaffolds}
        mv results/03_Assembly/{wildcards.sample}_assembly_dir/assembly_graph.fastg {output.graph}
        mv results/03_Assembly/{wildcards.sample}_assembly_dir/spades.log {output.spades_log}
        # rm -rf results/03_Assembly/{wildcards.sample}_assembly_dir
        """

rule run_quast:
    input:
        scaffolds_unfiltered = "results/03_Assembly/{sample}/{sample}_scaffolds.fasta",
        scaffolds = "results/03_Assembly/{sample}/{sample}_scaffolds_1k_filtered.fasta",
    output:
        quast_marker = "results/03_Assembly/{sample}/{sample}_quast_marker.done",
    params:
        quast_report = "results/03_Assembly/{sample}/{sample}_quast_report",
        jobname="{sample}_quast",
        account="pengel_spirit",
        runtime_s=convertToSec("0-07:10:00"),
    resources:
        mem_mb = convertToMb("50G")
    threads: 4
    log: "results/03_Assembly/{sample}/{sample}_quast.log"
    benchmark: "results/03_Assembly/{sample}/{sample}_quast.benchmark"
    conda: "../config/envs/trim-qc-env.yaml"
    shell:
        """
        quast.py -t {threads} -o {params.quast_report} {input.scaffolds} {input.scaffolds_unfiltered} &> {log}
        touch {output.quast_marker}
        """

rule filter_rename_scaffolds:
    input:
        scaffolds = "results/03_Assembly/{sample}/{sample}_scaffolds.fasta",
    output:
        scaffolds_filtered = "results/03_Assembly/{sample}/{sample}_scaffolds_1k_filtered.fasta",
    params:
        jobname="{sample}_filter_scaffold",
        account="pengel_spirit",
        runtime_s=convertToSec("0-00:10:00"),
        length_threshold = 1000,
        cov_threshold = 1,
    resources:
        mem_mb = convertToMb("50G")
    threads: 2
    log: "results/03_Assembly/{sample}/{sample}_filter_scaffold.log"
    benchmark: "results/03_Assembly/{sample}/{sample}_filter_scaffold.benchmark"
    conda: "../config/envs/mapping-bowtie-env.yaml"
    shell:
        """
        python3 scripts/filter_rename_scaffolds.py --in {input.scaffolds} --out {output.scaffolds_filtered} \
            --length_threshold {params.length_threshold} --cov_threshold {params.cov_threshold} \
            --sample {wildcards.sample} > {log}
        """