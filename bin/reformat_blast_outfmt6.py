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

if sys.argv[1] == '-v':
    print('1.0.0')
else:
    in_path = sys.argv[1]

in_data = gpf.ll(in_path)

for line in in_data:
    split_line = line.split()
    assert len(split_line) == 14
    output_line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(split_line[0], split_line[4], split_line[5], split_line[6], split_line[7], split_line[8], split_line[9], split_line[10], split_line[11], split_line[12], split_line[13], split_line[2])
    print(output_line)
