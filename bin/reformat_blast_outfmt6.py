#!/usr/bin/env python3
"""
Script for converting BLAST output from the format used by
BlobToolKit to standard outfmt 6

Written by Eerik Aunin @eeaunin

Adapted by Damon-Lee Pointon @DLBPointon
"""

import general_purpose_functions as gpf
import sys

# input format: qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue
# output format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

if sys.argv[1] == "-v":
    print("1.0.0")
    sys.exit()
else:
    in_path = sys.argv[1]

in_data = gpf.ll(in_path)

for line in in_data:
    s = line.split()
    assert len(s) == 14
    print(f"{s[0]}\t{s[4]}\t{s[5]}\t{s[6]}\t{s[7]}\t{s[8]}\t{s[9]}\t{s[10]}\t{s[11]}\t{s[12]}\t{s[13]}\t{s[2]}")
