"""
Code: Fragment Size Extraction (Version 3.2)
By: Songphon Sutthitthasakul (Moon)
Date: 30-05-2025
"""

import __init__ 
import os
import numpy as np
import pandas as pd
from subprocess import run

class cfMex_FS:
    def __init__(self, Input_path):
        self.max_size = 400
        self.Input_path = Input_path
        self.fragsize = np.zeros(self.max_size)

    def count_reads(self, mapq=30):
        df = pd.read_csv(
            self.Input_path,
            sep="\t",
            compression="gzip",
            header=None,
            names=["chr", "start", "end", "mapq", "strand"]
        )
        df_clean = df.loc[df["mapq"] >= mapq]
        df_clean["length"] = df_clean["end"] - df_clean["start"]
        df_clean = df_clean.loc[df_clean['length'] <= self.max_size]

        count_tmp = df_clean["length"].value_counts()
        for i in count_tmp.index:
            self.fragsize[i-1] =  count_tmp[i]

    def save_csv(self, Id, Output):
        with open(f'{Output}.csv', 'a') as outfile:
            run(['echo', f"{Id},{','.join(map(str,self.fragsize.tolist()))}"], stdout=outfile)

def FS_Run(Input, Id, Output):
    new_FS = cfMex_FS(Input)
    new_FS.count_reads()
    new_FS.save_csv(Id, Output)
