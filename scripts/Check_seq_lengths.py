import os
import glob
import gzip
from itertools import chain
from Bio import SeqIO
import numpy as np

SAMPLES = [os.path.basename(x) for x in glob.glob('results/00_RawData/from_Novogene/*') if os.path.basename(x).startswith('DNA')]
raw_paths_dict = {x: glob.glob(f'results/00_RawData/from_Novogene/{x}/*.fq.gz') for x in SAMPLES}
raw_files = list(chain(*raw_paths_dict.values()))

def remove_run_name_fluff(string):
    '''
    This function removes the run name fluff from the file name
    '''
    start = string.find('EKDN')
    end = string.find('SXC_')
    return(string[:start] + string[end+4:])

# median_lengths = {file: 0 for file in raw_files}
# for file in raw_files:
#     lengths = []
#     with gzip.open(file, "rt") as handle:
#         for record in SeqIO.parse(handle, "fastq"):
#             lengths.append(len(record.seq))
#             # print(len(record.seq))
#     print(f'{file}: {np.median(lengths)}')
#     median_lengths[file] = np.median(lengths)

# trimmed_paths_dict = {x: [os.path.join(f'results/01_TrimmingFiltering', remove_run_name_fluff(os.path.basename(y))) for y in glob.glob(f'results/00_RawData/from_Novogene/{x}/*.fq.gz')] for x in SAMPLES}
# trimmed_files = list(chain(*trimmed_paths_dict.values()))

# trimmed_median_lengths = {file: 0 for file in trimmed_files}
# for file in trimmed_files:
#     lengths = []
#     with gzip.open(file, "rt") as handle:
#         for record in SeqIO.parse(handle, "fastq"):
#             lengths.append(len(record.seq))
#             # print(len(record.seq))
#     trimmed_median_lengths[file] = np.median(lengths)

# with open('results/visualization/seq_lengths.txt', 'w') as f:
#     f.write('file,median_length\n')
#     for file, length in median_lengths.items():
#         f.write(f'{file},{length}\n')

# with open('results/visualization/trimmed_seq_lengths.txt', 'w') as f:
#     f.write('file,median_length\n')
#     for file, length in trimmed_median_lengths.items():
#         f.write(f'{file},{length}\n')

# read number of reads from fastq text files and make a table for R which is written in results/visualization
'''
There are multiple lanes for some samples and R1 and R2 reads for each lane.
The raw and trimmed read numbers are taken from the fastqc files made for each lane of each read set of each sample.
For host mapping all lanes are taken together by bowtie2 for each sample and the number mapped includes alignments from forward and reverse reads.
For host filtering done by bbmap it takes in all the lanes at once and reports mapping forwards and reverse reads separately.
'''
with open('results/visualization/number_of_reads.tsv', 'w+') as f_out:
    success = f_out.write(f'file\tsample\traw_reads\ttrimmed_reads\n')
    string_write = ''
    raw_reads = 0
    trimmed_reads = 0
    for file_path in raw_files:
        file = remove_run_name_fluff(os.path.basename(file_path))
        sample = file.split('_L')[0]
        with open(f'results/00_RawQC/{file.replace(".fq.gz", "_fastqc.txt")}', 'r') as f:
            for line in f:
                if line.startswith('Total Sequences'):
                    raw_reads = int(line.split('\t')[1].strip())
        try:
            with open(f'results/01_TrimmedQC/{file.replace(".fq.gz", "_fastqc.txt")}', 'r') as f:
                for line in f:
                    if line.startswith('Total Sequences'):
                        trimmed_reads = int(line.split('\t')[1].strip())
        except FileNotFoundError:
            trimmed_reads = 0
            print(f'No trimmed file for {file}')
        string_write = f'{file}\t{sample}\t{raw_reads}\t{trimmed_reads}\n'
        success = f_out.write(string_write)

with open('results/visualization/number_of_reads_mapping.tsv', 'w+') as f_out:
    success = f_out.write(f'sample\ttotal_reads\thost_reads\tpercent_host_relaxed\thost_mapped_depletion\tpercent_host_strict\tfinal_clean_reads\tpercent_final\n')
    for sample in SAMPLES:
        string_write = ''
        host_reads = 0
        total_reads = 0
        host_mapped_depletion = 0
        final_reads = 0
        percent_host_relaxed = 0
        percent_host_strict = 0
        percent_final = 0
        try:
            with open(f'results/02_HostMapping/{sample}/{sample}_bowtie_flagstat.tsv') as f:
                for line in f:
                    if 'total' in line:
                        total_reads = int(line.split()[0])
                    if 'primary mapped' in line:
                        host_reads = int(line.split()[0])
        except FileNotFoundError:
            host_reads = 0
            print(f'No host file for {sample}')
        try:
            with open(f'results/02_HostDepleting/{sample}/{sample}_host_depleting.log') as f:
                for line in f:
                    if 'mapped' in line:
                        host_mapped_depletion += int(line.split()[2])
        except FileNotFoundError:
            host_mapped_depletion = 0
            print(f'No host file for {sample}')
        final_reads = total_reads - host_mapped_depletion
        percent_host_relaxed = round(host_reads / total_reads * 100, 2)
        percent_host_strict = round(host_mapped_depletion / total_reads * 100, 2)
        percent_final = round(final_reads / total_reads * 100, 2)
        string_write = f'{sample}\t{total_reads}\t{host_reads}\t{percent_host_relaxed}\t{host_mapped_depletion}\t{percent_host_strict}\t{final_reads}\t{percent_final}\n'
        success = f_out.write(string_write)