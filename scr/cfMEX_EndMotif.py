"""
Code: End Motif Extraction (Version 2.3)
By: Songphon Sutthitthasakul (Moon)
Date: 29-05-2025
"""

import __init__ 
import gzip
import pysam
import itertools
import pandas as pd
from Bio import SeqIO
from subprocess import run
from collections import Counter
from .Utility import reverse_complement

class cfMex_EM:
    def __init__(self):
        self.reference = SeqIO.to_dict(SeqIO.parse(__init__.Reference_Genome, 'fasta'))
        self.nuc_dict = {}
    
    def count_reads(self, loc, Input_path, Id):
        nuc_comb = [''.join(p) for p in itertools.product(['A','C','G','T'], repeat=4)]

        bedfile = pysam.TabixFile(Input_path)
        EM_list =[]

        with gzip.open(Input_path, 'rt') as file:
            first_line = file.readline()
            if 'chr' not in first_line:
                loc[0] = loc[0].replace('chr','')

        for record in bedfile.fetch(loc[0], loc[1], loc[2]):
            fields = record.strip().split('\t')
            chrom, start, end, mapq, strand = fields[0], int(fields[1]), int(fields[2]), int(fields[3]), fields[4]

            if 'chr' not in chrom:
                chrom = f'chr{chrom}'

            if strand == '+' or strand == '0':
                record = self.reference[chrom]
                EM_list.append(str(record.seq[start:start+4]))
            elif strand == '-' or strand == '1':
                record = self.reference[chrom]
                EM_list.append(reverse_complement(record.seq[end-4:end]))

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


def EM_Run(loc, Input, Id):
    for locate in loc:
        new_EM = cfMex_EM()
        new_EM.count_reads(locate, Input, Id)

