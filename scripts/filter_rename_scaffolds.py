#!/usr/bin/env python3

import sys
import os
import argparse
import gzip

'''
usage:
    scripts/filter_rename_scaffolds.py --in {input.scaffolds} --out {output.scaffolds_filtered} &> {log}
'''

args = sys.argv
parser = argparse.ArgumentParser(description='Filter and rename scaffolds by appending sample name at the beginning of the header. Only fasta files from spades output.')
parser.add_argument('--input', type=str, help='Input scaffolds file', action='store')
parser.add_argument('--output', type=str, help='Output filtered scaffolds file', action='store')
parser.add_argument('--sample', type=str, help='Sample name to append to the header', action='store')
parser.add_argument('--length_threshold', type=int, help='Minimum length of contig to keep', action='store')
parser.add_argument('--cov_threshold', type=float, help='Minimum coverage of contig to keep', action='store')
args = parser.parse_args()

infile = args.input
outfile = args.output
sample = args.sample
length_threshold = args.length_threshold
cov_threshold = args.cov_threshold

def accept_contig(header):
    length = header.split("_")[3]
    cov = header.split("_")[5]
    if float(length) > length_threshold and float(cov) > cov_threshold:
        return(True)
    else:
        return(False)

write_line = False
count = 0
count_filtered = 0

with open(infile, 'r') as f_in:
    with open(outfile, 'w') as f_out:
        print("Filtering contigs from " + infile + " to " + outfile)
        for line in f_in:
            if line.startswith('>'):
                count += 1
                if accept_contig(line):
                    write_line = True
                    line = f'>{sample}--{line[1:]}'
                    count_filtered += 1
                else:
                    write_line = False
            if write_line:
                f_out.write(line)

print(f'Total contigs: {count}\n Contigs filtered: {count_filtered}')

    
