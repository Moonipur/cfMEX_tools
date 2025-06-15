"""
Code: Each Range of Fragment Size Extraction (Version 1.0)
By: Songphon Sutthitthasakul (Moon)
Date: 09-06-2025
"""

import os
from subprocess import run

import numpy as np
import pandas as pd

import __init__


class cfMex_Range:
    def __init__(self, Input_path, Format):
        self.Input_path = Input_path
        self.format = Format

    def split_fragment(self, loc, mapq=30):
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
        df_clean = df_clean.loc[df_clean["length"] > loc[0]]
        df_clean = df_clean.loc[df_clean["length"] <= loc[1]]

        df_clean = df_clean[["chr", "start", "end", "mapq", "strand"]]

        output = self.Input_path.replace(".hg38.bed.gz", f".{loc[0]}_{loc[1]}.hg38.bed")
        df_clean.to_csv(
            output,
            index=False,
            header=False,
            sep="\t",
        )


def Range_Run(loc, Input, Format):
    for cfRange in loc:
        new_FS_EM = cfMex_Range(Input, Format)
        new_FS_EM.split_fragment(cfRange)
