import os
import glob
import gzip
from itertools import chain
from Bio import SeqIO

raw_paths_dict = {x: glob.glob(f'results/00_RawData/from_Novogene/{x}/*.fq.gz') for x in SAMPLES}
raw_files = list(chain(*raw_paths_dict.values()))

median_lengths = {file: 0 for file in raw_files}
for file in raw_files:
    lengths = []
    with gzip.open(file, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            lengths.append(len(record.seq))
            print(len(record.seq))
    print(f'{file}: {np.medians(lengths)}')
    median_lengths[file] = np.median(lengths)

trimmed_paths_dict = {x: [os.path.join(f'results/01_TrimmingFiltering', remove_run_name_fluff(os.path.basename(y))) for y in glob.glob(f'results/00_RawData/from_Novogene/{x}/*.fq.gz')] for x in SAMPLES}
trimmed_files = list(chain(*trimmed_paths_dict.values()))

trimmed_median_lengths = {file: 0 for file in trimmed_files}
for file in trimmed_files:
    lengths = []
    with gzip.open(file, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            lengths.append(len(record.seq))
    trimmed_median_lengths[file] = np.median(lengths)

