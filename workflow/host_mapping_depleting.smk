#!/usr/bin/env python

"""
name: host_mapping_depleting
description: 
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
    log: "results/02_HostMapping/bowtie_index.log"
    benchmark: "results/02_HostMapping/bowtie_index.benchmark"
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
        bam = temp("results/02_HostMapping/{sample}/{sample}_bowtie.bam"),
        flagstat = "results/02_HostMapping/{sample}/{sample}_bowtie_flagstat.tsv",
    params:
        jobname="{sample}_map_to_host_bowtie",
        account="pengel_spirit",
        runtime_s=convertToSec("0-17:10:00"),
        input_files_string1 = lambda wildcards: ','.join(get_trimmed_files(sample = wildcards.sample, read = "1")),
        input_files_string2 = lambda wildcards: ','.join(get_trimmed_files(sample = wildcards.sample, read = "2")),
    resources:
        mem_mb = convertToMb("50G")
    threads: 4
    log: "results/02_HostMapping/{sample}/{sample}_host_mapping_bowtie.log"
    benchmark: "results/02_HostMapping/{sample}/{sample}_host_mapping_bowtie.benchmark"
    conda: "../config/envs/mapping-bowtie-env.yaml"
    shell:
        """
        bowtie2 -x {input.host_database} -1 {params.input_files_string1} \
            -2 {params.input_files_string1} \
            | samtools view -bh - | samtools sort - > {output.bam}
        samtools flagstat {output.bam} > {output.flagstat}
        """

# 02_HostMapping/{sample}_R1_filt.fq.gz and 02_HostMapping/{sample}_R2_filt.fq.gz are all reads not host-removed reads
# but 02_HostDepleting/{sample}_R1_filt.fq.gz and 02_HostDepleting/{sample}_R2_filt.fq.gz are host-depleted reads
# not trying 

# gets all host reads but might include some bacterial reads so too srtict for bacterial but very forgiving for host
rule get_non_host_reads:
    input:
        bam = "results/02_HostMapping/{sample}/{sample}_bowtie.bam",
    output:
        reads1 = "results/02_HostMapping/{sample}/{sample}_R1_filt.fq.gz", # all reads
        reads2 = "results/02_HostMapping/{sample}/{sample}_R2_filt.fq.gz", # all reads
        reads1_host = "results/02_HostMapping/{sample}/{sample}_R1_hostreads.fq.gz", # all host-mapped reads
        reads2_host = "results/02_HostMapping/{sample}/{sample}_R2_hostreads.fq.gz", # all host-mapped reads
        bam_unmapped = temp("results/02_HostMapping/{sample}/{sample}_host_unmapped.bam"), # not really used
        bam_host = temp("results/02_HostMapping/{sample}/{sample}_host.bam"),
    params:
        jobname="{sample}_get_non_host_reads",
        account="pengel_spirit",
        runtime_s=convertToSec("0-07:10:00"),
    resources:
        mem_mb = convertToMb("50G")
    threads: 4
    log: "results/02_HostMapping/{sample}/{sample}_get_non_host_reads.log"
    benchmark: "results/02_HostMapping/{sample}/{sample}_get_non_host_reads.benchmark"
    conda: "../config/envs/mapping-bowtie-env.yaml"
    shell:
        """
        samtools view -b -f 4 {input.bam} | samtools sort -n - > {output.bam_unmapped}
        samtools view -b -F 4 {input.bam} | samtools sort -n - > {output.bam_host}
        outreads1={output.reads1}
        outreads2={output.reads2}
        bedtools bamtofastq -i {input.bam} -fq ${{outreads1/.gz/}} -fq2 ${{outreads2/.gz/}} &> {log}
        gzip ${{outreads1/.gz/}}
        gzip ${{outreads2/.gz/}}
        outreads_h_1={output.reads1_host}
        outreads_h_2={output.reads2_host}
        bedtools bamtofastq -i {output.bam_host} -fq ${{outreads_h_1/.gz/}} -fq2 ${{outreads_h_2/.gz/}} &> {log}
        gzip ${{outreads_h_1/.gz/}}
        gzip ${{outreads_h_2/.gz/}}
        """

rule bbmap_index_host:
    input:
        host_genome = "data/host_genome/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
    output:
        bbmap_index_done = "data/host_genome/bbmap_index.done"
    params:
        bbmap_index = "data/host_genome/",
        jobname="bbmap_index_host",
        account="pengel_spirit",
        runtime_s=convertToSec("0-07:10:00"),
    resources:
        mem_mb = convertToMb("50G")
    threads: 4
    log: "results/02_HostDepleting/bbmap_index_host.log"
    benchmark: "results/02_HostDepleting/bbmap_index_host.benchmark"
    conda: "../config/envs/mapping-bowtie-env.yaml"
    shell:
        """
        bbmap.sh -Xmx25g ref={input.host_genome} path={params.bbmap_index} &> {log}
        touch {output.bbmap_index_done}
        """

# gets all bacterial reads but might include some host reads so too srtict for host but very forgiving for bacterial will be used for further analysis
rule host_depleting_bbmap:
    input:
        reads1 = lambda wildcards: get_trimmed_files(sample = wildcards.sample, read = "1"),
        reads2 = lambda wildcards: get_trimmed_files(sample = wildcards.sample, read = "2"),
        bbmap_index = "data/host_genome/bbmap_index.done"
    output:
        reads1 = "results/02_HostDepleting/{sample}/{sample}_R1_filt.fq.gz",
        reads2 = "results/02_HostDepleting/{sample}/{sample}_R2_filt.fq.gz",
    params:
        host_genome_path = "data/host_genome",
        input_files_string1 = lambda wildcards: order_csv_list_string(','.join(get_trimmed_files(sample = wildcards.sample, read = "1"))),
        input_files_string2 = lambda wildcards: order_csv_list_string(','.join(get_trimmed_files(sample = wildcards.sample, read = "2"))),
        jobname="{sample}_host_depleting",
        account="pengel_spirit",
        runtime_s=convertToSec("0-15:10:00"),
    resources:
        mem_mb = convertToMb("50G")
    threads: 4
    log: "results/02_HostDepleting/{sample}/{sample}_host_depleting.log"
    benchmark: "results/02_HostDepleting/{sample}/{sample}_host_depleting.benchmark"
    conda: "../config/envs/mapping-bowtie-env.yaml"
    shell:
        """
        bbwrap.sh -Xmx25g threads={threads} minid=0.95 maxindel=3 \
            bwr=0.16 bw=12 quickmatch fast minhits=2 path={params.host_genome_path} \
            in1={params.input_files_string1} in2={params.input_files_string2} outu1={output.reads1} outu2={output.reads2} 2>> {log}
        """

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