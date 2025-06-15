"""
Code: BAM-2-Fragment Convertion (Version 2.2)
By: Songphon Sutthitthasakul (Moon)
Date: 12-06-2025
"""

import itertools
import sys
from collections import Counter

import numpy as np
import pandas as pd

# from .Utility import reverse_complement


class BAM2FRAG:
    def __init__(self, input, ref_ver="hg38"):
        self.Input = input
        self.mapq = 30
        self.ref = ref_ver

    def Convert(self):
        dummy_chrom_hg19 = {
            "1": 1,
            "2": 2,
            "3": 3,
            "4": 4,
            "5": 5,
            "6": 6,
            "7": 7,
            "8": 8,
            "9": 9,
            "10": 10,
            "11": 11,
            "12": 12,
            "13": 13,
            "14": 14,
            "15": 15,
            "16": 16,
            "15": 15,
            "16": 16,
            "17": 17,
            "18": 18,
            "19": 19,
            "20": 20,
            "21": 21,
            "22": 22,
            "X": 23,
            "Y": 24,
        }

        dummy_chrom_hg38 = {
            "chr1": 1,
            "chr2": 2,
            "chr3": 3,
            "chr4": 4,
            "chr5": 5,
            "chr6": 6,
            "chr7": 7,
            "chr8": 8,
            "chr9": 9,
            "chr10": 10,
            "chr11": 11,
            "chr12": 12,
            "chr13": 13,
            "chr14": 14,
            "chr15": 15,
            "chr16": 16,
            "chr15": 15,
            "chr16": 16,
            "chr17": 17,
            "chr18": 18,
            "chr19": 19,
            "chr20": 20,
            "chr21": 21,
            "chr22": 22,
            "chrX": 23,
            "chrY": 24,
        }

        data = pd.read_csv(
            self.Input,
            sep="\t",
            compression="gzip",
            header=None,
            names=["chr", "start1", "end1", "start2", "end2", "mapq", "strand"],
        )

        data = data.loc[data["mapq"] >= self.mapq]

        direct = data["strand"].values.tolist()
        start1 = data["start1"].values.tolist()
        start2 = data["start2"].values.tolist()
        end1 = data["end1"].values.tolist()
        end2 = data["end2"].values.tolist()

        data = data[["chr", "mapq"]]

        start = []
        end = []
        strand = []

        for idx, i in enumerate(direct):
            start.append(np.min([start1[idx], start2[idx]]))
            end.append(np.max([end1[idx], end2[idx]]))
            if i == "+":
                strand.append(0)
            else:
                strand.append(1)

        data["start"] = start
        data["end"] = end
        data["strand"] = strand

        data = data[["chr", "start", "end", "mapq", "strand"]]

        if self.ref == "hg38":
            data["chr"] = dummy_chrom_hg38[f"{data.iloc[0, 0]}"]
        elif self.ref == "hg19":
            data["chr"] = dummy_chrom_hg19[f"{data.iloc[0, 0]}"]

        data = data.sort_values(by=["start", "end"])

        output = self.Input.replace(".frag.gz", f".{self.ref}.bed")
        data.to_csv(
            output,
            index=False,
            header=False,
            sep="\t",
        )


if __name__ == "__main__":
    BAM2FRAG(sys.argv[1], sys.argv[2]).Convert()
