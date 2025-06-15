"""
Code: End Motif Extraction (Version 2.4)
By: Songphon Sutthitthasakul (Moon)
Date: 06-06-2025
"""

import gzip
import itertools
from collections import Counter
from subprocess import run

import pandas as pd
import pysam
from Bio import SeqIO

import __init__

from .Utility import reverse_complement


class cfMex_EM:
    def __init__(self, Format="c5"):
        self.reference = SeqIO.to_dict(SeqIO.parse(__init__.Reference_Genome, "fasta"))
        self.nuc_dict_ref = {}
        self.nuc_dict_ori = {}
        self.format = Format

    def count_reads(self, loc, Input_path, Id):
        nuc_comb = [
            "".join(p) for p in itertools.product(["A", "C", "G", "T"], repeat=4)
        ]

        bedfile = pysam.TabixFile(Input_path)
        EM_ref = []
        EM_ori = []

        with gzip.open(Input_path, "rt") as file:
            first_line = file.readline()
            if "chr" not in first_line:
                loc[1] = loc[1].replace("chr", "")

        for record in bedfile.fetch(loc[1], loc[2], loc[3]):
            fields = record.strip().split("\t")
            if self.format == "c5":
                chrom, start, end, mapq, strand = (
                    fields[0],
                    int(fields[1]),
                    int(fields[2]),
                    int(fields[3]),
                    fields[4],
                )
            elif self.format == "c6":
                chrom, start, end, mapq, motif, strand = (
                    fields[0],
                    int(fields[1]),
                    int(fields[2]),
                    int(fields[3]),
                    int(fields[4]),
                    fields[5],
                )
                EM_ori.append(motif)

            if "chr" not in chrom:
                chrom = f"chr{chrom}"

            if strand == "+" or strand == "0":
                record = self.reference[chrom]
                EM_ref.append(str(record.seq[start : start + 4]))
            elif strand == "-" or strand == "1":
                record = self.reference[chrom]
                EM_ref.append(reverse_complement(record.seq[end - 4 : end]))

        count_ref = Counter(EM_ref)
        if self.format == "c6":
            count_ori = Counter(EM_ori)

        for nix, nuc in enumerate(nuc_comb):
            self.nuc_dict_ref[nuc] = count_ref[nuc]
            if self.format == "c6":
                self.nuc_dict_ori[nuc] = count_ori[nix + 1]

        data_ref = map(str, list(self.nuc_dict_ref.values()))
        if self.format == "c6":
            data_ori = map(str, list(self.nuc_dict_ori.values()))

        with open(f"{Id}_Metadata_ref.csv", "a") as meta:
            run(["echo", f'{",".join(data_ref)}'], stdout=meta)

        if self.format == "c6":
            with open(f"{Id}_Metadata_ori.csv", "a") as meta:
                run(["echo", f'{",".join(data_ori)}'], stdout=meta)

    @staticmethod
    def save_csv(Id, Output, Format):
        data = pd.read_csv(f"{Id}_Metadata_ref.csv", header=None)
        sum_data = map(str, data.sum(axis=0).values.tolist())
        with open(f"{Output}.csv", "a") as outfile:
            run(["echo", f'{Id},{",".join(sum_data)}'], stdout=outfile)

        if Format == "c6":
            data = pd.read_csv(f"{Id}_Metadata_ori.csv", header=None)
            sum_data = map(str, data.sum(axis=0).values.tolist())
            with open(f"{Output}_ori.csv", "a") as outfile:
                run(["echo", f'{Id},{",".join(sum_data)}'], stdout=outfile)


def EM_Run(loc, Input, Id, Format):
    for locate in loc:
        new_EM = cfMex_EM(Format)
        new_EM.count_reads(locate, Input, Id)
