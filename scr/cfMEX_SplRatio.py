"""
Code: Short per Long Ratio Extraction (Version 1.3)
By: Songphon Sutthitthasakul (Moon)
Date: 12-06-2025
"""

import gzip
from subprocess import run

import pandas as pd
import pysam

import __init__


class cfMex_SR:
    def __init__(self, Format="c5"):
        self.nuc_dict = {}
        self.format = Format

    def count_reads(self, loc, Input_path, Id):
        short = 0
        long = 0

        bedfile = pysam.TabixFile(Input_path)

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

            if "chr" not in chrom:
                chrom = f"chr{chrom}"

            if end - start >= 151 and end - start <= 250 and mapq >= 30:
                long += 1
            elif end - start >= 90 and end - start <= 150 and mapq >= 30:
                short += 1

        if long == 0:
            spl = 0
        else:
            spl = short / long

        with open(f"{Id}_Metadata.csv", "a") as meta:
            run(["echo", f"{loc[0]},{spl}"], stdout=meta)

    @staticmethod
    def save_csv(Id, Output, loc_n):
        data = pd.read_csv(f"{Id}_Metadata.csv", header=None, names=["index", "spl"])
        data = data.sort_values(by=["index"])
        if data.shape[0] == loc_n:
            sum_data = ",".join(map(str, data["spl"].tolist()))
            with open(f"{Output}.csv", "a") as outfile:
                run(["echo", f"{Id},{sum_data}"], stdout=outfile)
        else:
            print(f"Numbers of Locations is wrong: {data.shape[0]}")


def SR_Run(loc, Input, Id, Format):
    for locate in loc:
        new_SR = cfMex_SR(Format)
        new_SR.count_reads(locate, Input, Id)
