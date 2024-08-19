Work in Progress! This pipeline is not intended to be run independently without prior set-up. If you would like to use it or its parts and have questions or require clarification, contact Aiswarya Prasad. (aiswarya.prasad@unil.ch)

# Comparison of wild and managed Apis mellifera

The aim of this pipeline is to document the steps used to process raw shotgun metagenomic data (R1 and R2 fastq reads), assemble and bin scaffolds into MAGs, cluster them into magOTUs, estimate the coverage of each of the magOTUs and then finally also SNP profiling across samples based on the high quality MAGs chosen. Downstream analysis is performed using independent scripts and documented in their respective directories. All scripts are present in the scripts directory in general.

# Directory structure

* config
* data
* results
* scripts
* workflow

# Other details

snakemake: 7.28.1

running using the command

```bash
snakemake -p -r \
--use-conda --conda-prefix /work/FAC/FBM/DMF/pengel/spirit/aprasad/snakemake-conda-envs-2024 \
--conda-frontend mamba --profile slurm --restart-times 0 \
--keep-going --rerun-incomplete --keep-going --rerun-incomplete \
--rerun-triggers mtime --cluster-cancel scancel --jobs 50 -n
```
