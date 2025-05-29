"""
Code: BAM-2-Fragment Convertion (Version 1.3)
By: Songphon Sutthitthasakul (Moon)
Date: 29-05-2025
"""

import numpy as np
import pandas as pd

class BAM2FRAG:
    def __init__(self, input):
        self.Input = input
        self.mapq = 30
        
    def Convert(self):
        dummy_chrom = {
            'chr1': 1, 'chr2': 2,
            'chr3': 3, 'chr4': 4,
            'chr5': 5, 'chr6': 6,
            'chr7': 7, 'chr8': 8,
            'chr9': 9, 'chr10': 10,
            'chr11': 11, 'chr12': 12,
            'chr13': 13, 'chr14': 14,
            'chr15': 15, 'chr16': 16,
            'chr15': 15, 'chr16': 16,
            'chr17': 17, 'chr18': 18,
            'chr19': 19, 'chr20': 20,
            'chr21': 21, 'chr22': 22,
            'chrX': 23, 'chrY': 24,
        } 

        data = pd.read_csv(
            self.Input,
            sep="\t",
            compression="gzip",
            header=None,
            names=["chr", "start", "mapq", "length"]
        )

        length = data['length'].values.tolist()
        start = data['start'].values.tolist()

        end = []
        strand = []

        for idx, i in enumerate(length):
            end.append(start[idx] + np.abs(length[idx]))
            if i >= 0:
                strand.append(0)
            else:
                strand.append(1)
        
        data['end'] = end
        data['strand'] = strand

        data = data[['chr', 'start', 'end', 'mapq', 'strand']]
        data = data.loc[data["mapq"] >= self.mapq]
        data['chr'] = dummy_chrom[data.iloc[0,0]]

        data = data.sort_values(by=['start', 'end'])

        output = self.Input.replace('.frag.gz','.hg38.bed')
        data.to_csv(
            output, 
            index=False, 
            header=False, 
            sep='\t',
        )
            
