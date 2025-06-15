"""
Code: BAM-2-Fragment Convertion (Version 1.4)
By: Songphon Sutthitthasakul (Moon)
Date: 06-06-2025
"""

import itertools
from collections import Counter

import numpy as np
import pandas as pd

from .Utility import reverse_complement


class BAM2FRAG:
    def __init__(self, input):
        self.Input = input
        self.mapq = 30

    def Convert(self):
        dummy_chrom = {
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
            names=["chr", "start", "mapq", "length", "motif"],
        )

        length = data["length"].values.tolist()
        start = data["start"].values.tolist()
        motif = data["motif"].values.tolist()

        nuc_comb = [
            "".join(p) for p in itertools.product(["A", "C", "G", "T"], repeat=4)
        ]
        nuc_dict = {}

        for nix, nuc in enumerate(nuc_comb):
            nuc_dict[nuc] = nix + 1

        end = []
        strand = []
        new_motif = []

        for idx, i in enumerate(length):
            end.append(start[idx] + np.abs(length[idx]))
            if i >= 0:
                strand.append(0)
                try:
                    new_motif.append(nuc_dict[f"{motif[idx][:4]}"])
                except:
                    new_motif.append(0)
            else:
                strand.append(1)
                try:
                    new_motif.append(nuc_dict[f"{reverse_complement(motif[idx][-4:])}"])
                except:
                    new_motif.append(0)

        data["end"] = end
        data["strand"] = strand
        data["motif"] = new_motif

        data = data[["chr", "start", "end", "mapq", "motif", "strand"]]
        data = data.loc[data["mapq"] >= self.mapq]
        data["chr"] = dummy_chrom[data.iloc[0, 0]]

        data = data.sort_values(by=["start", "end"])

        output = self.Input.replace(".frag.gz", ".hg38.bed")
        data.to_csv(
            output,
            index=False,
            header=False,
            sep="\t",
        )
