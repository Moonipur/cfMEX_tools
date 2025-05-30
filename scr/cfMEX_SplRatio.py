"""
Code: Short per Long Ratio Extraction (Version 1.1)
By: Songphon Sutthitthasakul (Moon)
Date: 30-05-2025
"""

import __init__ 
import gzip
import pysam
import pandas as pd
from subprocess import run

class cfMex_SR:
    def __init__(self):
        self.nuc_dict = {}
    
    def count_reads(self, loc, Input_path, Id):
        short = 0
        long = 0

        bedfile = pysam.TabixFile(Input_path)

        with gzip.open(Input_path, 'rt') as file:
            first_line = file.readline()
            if 'chr' not in first_line:
                loc[1] = loc[1].replace('chr','')
                
        for record in bedfile.fetch(loc[1], loc[2], loc[3]):
            fields = record.strip().split('\t')
            chrom, start, end, mapq, strand = fields[0], int(fields[1]), int(fields[2]), int(fields[3]), fields[4]

            if 'chr' not in chrom:
                chrom = f'chr{chrom}'

            if end - start >= 151 and end - start <= 250:
                long += 1
            elif end - start >= 90 and end - start <= 150:
                short += 1

        if long == 0:
            spl = 0
        else:
            spl = short / long

        with open(f'{Id}_Metadata.csv','a') as meta:
            run(['echo', f'{loc[0]},{spl}'], stdout=meta)

    @staticmethod
    def save_csv(Id, Output):
        data = pd.read_csv(f'{Id}_Metadata.csv', header=None, names=['index','spl'])
        data = data.sort_values(by=['index'])
        sum_data = ','.join(map(str,data['spl'].tolist()))
        with open(f'{Output}.csv','a') as outfile:
            run(['echo', sum_data], stdout=outfile)


def SR_Run(loc, Input, Id):
    for locate in loc:
        new_SR = cfMex_SR()
        new_SR.count_reads(locate, Input, Id)
