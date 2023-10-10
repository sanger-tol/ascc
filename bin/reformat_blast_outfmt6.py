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
    print('\t'.join( s[0:1] + s[4:] + s[2:3] ))
