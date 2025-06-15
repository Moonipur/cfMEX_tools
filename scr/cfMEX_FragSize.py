"""
Code: Fragment Size Extraction (Version 3.3)
By: Songphon Sutthitthasakul (Moon)
Date: 06-06-2025
"""

import os
from subprocess import run

import numpy as np
import pandas as pd

import __init__


class cfMex_FS:
    def __init__(self, Input_path, Format):
        self.max_size = 400
        self.Input_path = Input_path
        self.fragsize = np.zeros(self.max_size)
        self.format = Format

    def count_reads(self, mapq=30):
        if self.format == "c5":
            df = pd.read_csv(
                self.Input_path,
                sep="\t",
                compression="gzip",
                header=None,
                names=["chr", "start", "end", "mapq", "strand"],
            )
        elif self.format == "c6":
            df = pd.read_csv(
                self.Input_path,
                sep="\t",
                compression="gzip",
                header=None,
                names=["chr", "start", "end", "mapq", "motif", "strand"],
            )
        df_clean = df.loc[df["mapq"] >= mapq]
        df_clean["length"] = df_clean["end"] - df_clean["start"]
        df_clean = df_clean.loc[df_clean["length"] <= self.max_size]

        count_tmp = df_clean["length"].value_counts()
        for i in count_tmp.index:
            self.fragsize[i - 1] = count_tmp[i]

    def save_csv(self, Id, Output):
        with open(f"{Output}.csv", "a") as outfile:
            run(
                ["echo", f"{Id},{','.join(map(str,self.fragsize.tolist()))}"],
                stdout=outfile,
            )


def FS_Run(Input, Id, Output, Format):
    new_FS = cfMex_FS(Input, Format)
    new_FS.count_reads()
    new_FS.save_csv(Id, Output)
