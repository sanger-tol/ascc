#!/usr/bin/env python3
"""
Script for filtering BLAST results to keep only top hits for each sequence
Argument1: path to BLAST output (format 6)
Output: BLAST results filtered to keep only the hits with top scores for each sequence

Written by Eerik Aunin @eeaunin

Adapted by Damon-Lee Pointon @DLBPointon
"""

import sys
import pandas as pd


if sys.argv[1] == "-v":
    print("1.0.0")
    sys.exit()
else:
    in_path = sys.argv[1]

df = pd.read_csv(in_path, header=None, sep="\t")

contigs = list(set(df[0].tolist()))
contigs.sort()

for contig in contigs:
    df2 = df.loc[df[0] == contig]
    max_score = max(df2[11])
    df3 = df2.loc[df2[11] == max_score]
    df3 = df3.iloc[0]
    df3_len = len(df3)
    for counter, item in enumerate(df3):
        item = str(item)
        if counter == df3_len - 1:
            print(item)
        else:
            print(item + "\t", end="")
