"""
Code: End Motif Extraction (Version 2.0)
By: Songphon Sutthitthasakul (Moon)
Date: 20-05-2025
"""

from __init__ import *
import os
import pysam
import argparse
import itertools
import pandas as pd
from Bio import SeqIO
from time import time
from subprocess import run
from collections import Counter
from multiprocessing import Pool

#################################################################################################################
# SETTING

Reference_Genome = '/mnt/sas/ref/hg38/v0/Homo_sapiens_assembly38.fasta'
EM_window = 'EM_Window.csv'

#################################################################################################################

class cfMex_EM:
    def __init__(self):
        self.reference = SeqIO.to_dict(SeqIO.parse(Reference_Genome, 'fasta'))
        self.nuc_dict = {}
    
    @staticmethod
    def reverse_complement(seq):
        seq = seq[::-1]
        new_seq = []
        for i in seq:
            if i == 'A':
                new_seq.append('T')
            elif i == 'T':
                new_seq.append('A')
            elif i == 'C':
                new_seq.append('G')
            elif i == 'G':
                new_seq.append('C')
            else:
                new_seq.append('N')
        return ''.join(new_seq)
    
    def count_reads(self, loc, Input_path, Id):
        nuc_comb = [''.join(p) for p in itertools.product(['A','C','G','T'], repeat=4)]

        bedfile = pysam.TabixFile(Input_path)
        EM_list =[]

        for record in bedfile.fetch(loc[0], loc[1], loc[2]):
            fields = record.strip().split('\t')
            chrom, start, end, mapq, strand = fields[0], int(fields[1]), int(fields[2]), int(fields[3]), fields[4]

            if strand == '+':
                record = self.reference[chrom]
                EM_list.append(str(record.seq[start:start+4]))
            elif strand == '-':
                record = self.reference[chrom]
                EM_list.append(cfMex_EM.reverse_complement(record.seq[end-4:end]))

        count = Counter(EM_list)
        for nuc in nuc_comb:
            self.nuc_dict[nuc] = count[nuc]

        data = map(str,list(self.nuc_dict.values()))
        with open(f'{Id}_Metadata.csv','a') as meta:
            run(['echo', f'{",".join(data)}'], stdout=meta)

    @staticmethod
    def save_csv(Id, Output):
        data = pd.read_csv(f'{Id}_Metadata.csv', header=None)
        sum_data = map(str,data.sum(axis=0).values.tolist())
        with open(f'{Output}.csv','a') as outfile:
            run(['echo', f'{Id},{",".join(sum_data)}'], stdout=outfile)



def EM_Run(loc):
    for locate in loc:
        new_EM = cfMex_EM()
        new_EM.count_reads(locate, args.input, args.id)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract End-motif feature from cfDNA Fragment file.")
    parser.add_argument("--input", '-i', type=str, help="Fragment file path", required=True)
    parser.add_argument("--output", '-o', type=str, help="Output file name", required=False, default='EM_output')
    parser.add_argument("--id", type=str, help="Sample ID", required=False)
    parser.add_argument("--thread", '-t', type=int, help="Number of Threads (default: 2)", required=False, default=2)
    args = parser.parse_args()

    start_time = time()

    location = Threads_Location(EM_window, args.thread)

    with Pool(args.thread) as p:
        p.map(EM_Run, location)
    cfMex_EM().save_csv(args.id, args.output)
    os.remove(f'{args.id}_Metadata.csv')
    
    Time_Stamp(start_time)

    
