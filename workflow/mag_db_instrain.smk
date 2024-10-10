#!/usr/bin/env python

"""
name: mag_db_instrain
description: Takes binning results from metabat2 and summarises it then runs checkm, gtdbtk, drep on mags and creates a filtered mag database for instrain (non-redundant) and redundant for mapping to
"""

rule make_mag_rep_database:
    input:
        collect_mags_marker = "results/05_MAGs_collection/collect_mags.done",
        mag_metadata_summary = lambda wildcards: checkpoints.mag_metadata_summary.get().output.metadata,
        rep_mags = lambda wildcards: expand("results/05_MAGs_collection/MAGs/{mag}.fa", mag = get_rep_mags(checkpoints.mag_metadata_summary.get().output.metadata)),
    output:
        mag_rep_database = "results/06_instrain/prepare_mags_db/mag_rep_database.fa"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-02:10:00"),
    threads: 4
    log: "results/06_instrain/prepare_mags_db/mag_rep_database.log"
    benchmark: "results/06_instrain/prepare_mags_db/mag_rep_database.benchmark"
    conda: "../config/envs/scripts-env.yaml"
    shell:
        """
        python scripts/make_mag_rep_database.py \
                --collect_mags_marker {input.collect_mags_marker} \
                --mag_metadata_summary {input.mag_metadata_summary} \
                --mag_rep_database {output.mag_rep_database}
        """

rule get_genes_mag_rep_database:
    input:
        mag_rep_database = "results/06_instrain/prepare_mags_db/mag_rep_database.fa"
    output:
        instrain_genes_faa = "results/06_instrain/prepare_mags_db/mag_rep_database_genes.faa",
        instrain_genes_fna = "results/06_instrain/prepare_mags_db/mag_rep_database_genes.fna"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-02:10:00"),
    threads: 4
    log: "results/06_instrain/prepare_mags_db/get_genes_mag_rep_database.log"
    benchmark: "results/06_instrain/prepare_mags_db/get_genes_mag_rep_database.benchmark"
    conda: "../config/envs/genes-env.yaml"
    shell:
        """
        prodigal -i {input.mag_rep_database} -d {output.instrain_genes_fna} -a {output.instrain_genes_faa} -p meta &> {log}
        """

rule bowtie_index_mags_db:
    input:
        mag_rep_database = "results/06_instrain/prepare_mags_db/mag_rep_database.fa"
    output:
        bowtie_index = directory("results/06_instrain/prepare_mags_db/mag_rep_database_bowtie_index")
    params:
        prefix = "mag_rep_database",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="bowtie_index_mags_db:",
        account="pengel_spirit",
        runtime_s=convertToSec("0-07:10:00"),
    resources:
        mem_mb = convertToMb("50G")
    threads: 4
    log: "results/06_instrain/prepare_mags_db/bowtie_index_mags_db:.log"
    benchmark: "results/06_instrain/prepare_mags_db/bowtie_index_mags_db:.benchmark"
    conda: "../config/envs/mapping-bowtie-env.yaml"
    shell:
        """
        mkdir -p {output.bowtie_index}
        bowtie2-build -f --threads {threads} {input.mag_rep_database} {output.bowtie_index}/{params.prefix}
        """

rule map_to_rep_MAGs_bowtie2:
    input:
        reads1 = "results/02_HostDepleting/{sample}/{sample}_R1_filt.fq.gz",
        reads2 = "results/02_HostDepleting/{sample}/{sample}_R2_filt.fq.gz",
        bowtie_index = "results/06_instrain/prepare_mags_db/mag_rep_database_bowtie_index",
        mag_rep_database = "results/06_instrain/prepare_mags_db/mag_rep_database.fa"
    output:
        bam = temp(SCRATCH_PATH+"06_instrain/bowtie2_mapping_temp/{sample}_bowtie.bam"),
        flagstat = "results/06_instrain/bowtie2_mapping/{sample}/{sample}_bowtie_flagstat.tsv",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_map_to_rep_MAGs_bowtie",
        account="pengel_spirit",
        runtime_s=convertToSec("0-17:10:00"),
    resources:
        mem_mb = convertToMb("50G")
    threads: 4
    log: "results/06_instrain/bowtie2_mapping/{sample}/{sample}_map_to_rep_MAGs_bowtie.log"
    benchmark: "results/06_instrain/bowtie2_mapping/{sample}/{sample}_map_to_rep_MAGs_bowtie.benchmark"
    conda: "../config/envs/mapping-bowtie-env.yaml"
    shell:
        """
        bowtie2 -X 1000 -x {input.mag_rep_database} -1 {input.reads1} -2 {input.reads2} | samtools view -bh - | samtools sort - > {output.bam}
        samtools flagstat {output.bam} > {output.flagstat}
        """

# A .text file with two columns separated by tabs, 
# where the first column is the name of a scaffold 
# and the second column is the name of the bin / genome the scaffold belongs to.
rule make_scaffold_to_bin_file:
    input:
        rep_mags_db = "results/06_instrain/prepare_mags_db/mag_rep_database.fa",
        rep_mags = lambda wildcards: expand("results/05_MAGs_collection/MAGs/{mag}.fa", mag = get_rep_mags(checkpoints.mag_metadata_summary.get().output.metadata)),
    output:
        scaffold_to_bin_file = "results/06_instrain/prepare_mags_db/scaffold_to_bin_file.tsv"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="make_scaffold_to_bin_file",
        account="pengel_spirit",
        runtime_s=convertToSec("0-01:10:00"),
    resources:
        mem_mb = convertToMb("8G")
    threads: 4
    log: "results/06_instrain/prepare_mags_db/make_scaffold_to_bin_file.log"
    benchmark: "results/06_instrain/prepare_mags_db/make_scaffold_to_bin_file.benchmark"
    run:
        with open(output.scaffold_to_bin_file, "w") as f:
            for mag in input.rep_mags:
                with open(mag, "r") as m:
                    for line in m:
                        mag_name = os.path.basename(mag).split(".")[0]
                        if line.startswith(">"):
                            scaffold = line.strip().split(">")[1]
                            f.write(f"{scaffold}\t{mag_name}\n")


rule instrain_profile_db_mode:
    input:
        bam = SCRATCH_PATH+"06_instrain/bowtie2_mapping_temp/{sample}_bowtie.bam",
        mag_rep_database = "results/06_instrain/prepare_mags_db/mag_rep_database.fa",
        scaffold_to_bin_file = "results/06_instrain/prepare_mags_db/scaffold_to_bin_file.tsv",
        instrain_genes_file = "results/06_instrain/prepare_mags_db/mag_rep_database_genes.fna",
    output:
        marker = touch("results/06_instrain/02_instrain_profile_db_mode/{sample}_profile_db_mode.done/")
    params:
        outdir = "results/06_instrain/02_instrain_profile_db_mode/{sample}",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_instrain_profile_db_mode",
        account="pengel_spirit",
        runtime_s=convertToSec("0-22:10:00"),
    resources:
        mem_mb = convertToMb("250G")
    threads: 8
    conda: "../config/envs/instrain_env.yaml"
    log: "results/06_instrain/02_instrain_profile_db_mode/{sample}_instrain_profile_db_mode.log"
    benchmark: "results/06_instrain/02_instrain_profile_db_mode/{sample}_instrain_profile_db_mode.benchmark"
    shell:
        """
        inStrain profile {input.bam} {input.mag_rep_database} -o {params.outdir} \
                -p {threads} -g {input.instrain_genes_file} \
                --max_insert_relative 5 --database_mode \
                -s {input.scaffold_to_bin_file}
        touch {output.marker}
        """

rule instrain_profile_plot:
    input:
        marker = "results/06_instrain/02_instrain_profile_db_mode/{sample}_profile.done/"
    output:
        done = touch("results/06_instrain/04_instrain_plot_marker/{sample}_profile_plots.done")
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        jobname="{sample}_instrain_profile_plot",
        account="pengel_spirit",
        runtime_s=convertToSec("0-07:10:00"),
    resources:
        mem_mb = convertToMb("100G")
    threads: 16
    log: "results/06_instrain/04_instrain_plot_marker/{sample}_profile_plots.log"
    benchmark: "results/06_instrain/04_instrain_plot_marker/{sample}_profile_plots.benchmark"
    conda: "../config/envs/instrain_env.yaml"
    shell:
        """
        profile={input.marker}
        profile=${{profile/_profile.done/}}
        inStrain plot -i ${{profile}} -pl a -p {threads}
        touch {output.done}
        """

rule instrain_compare:
    input:
        markers = expand("results/06_instrain/02_instrain_profile_db_mode/{sample}_profile_db_mode.done/", sample=SAMPLES),
        scaffold_to_bin_file = "results/06_instrain/prepare_mags_db/scaffold_to_bin_file.tsv"
    output:
        compare_marker = touch("results/06_instrain/03_instrain_compare/all_compared.done"),
    params:
        outdir = lambda wildcards: f"results/06_instrain/03_instrain_compare/20230313_MAGs_rep_db",
        profiles = [f"results/06_instrain/02_instrain_profile_db_mode/{sample}" for sample in SAMPLES],
        mailto="aiswarya.prasad@unil.ch",
        account="pengel_spirit",
        runtime_s=convertToSec("3-00:00:00"),
    resources:
        mem_mb = convertToMb("200G")
    conda: "../config/envs/instrain_env.yaml"
    log: "results/06_instrain/03_instrain_compare/log_files/instrain_compare_20230313_MAGs_rep_db.log"
    benchmark: "results/06_instrain/03_instrain_compare/log_files/instrain_compare_20230313_MAGs_rep_db.benchmark"
    threads: 8
    shell:
        """
        inStrain compare -i {params.profiles} -s {input.scaffold_to_bin_file} \
                -p {threads} -o {params.outdir} \
                --database_mode || touch {output.compare_marker}.singleton
        touch {output.compare_marker}
        """