#!/usr/bin/env python3
"""
Script reformatting Diamond BLASTX results so that they can used to make a BlobToolKit dataset
input:
--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms sphylums salltitles
output:
--outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

Written by Eerik Aunin @eeaunin

Adapted by Damon-Lee Pointon @DLBPointon
"""

import sys
import general_purpose_functions as gpf

if sys.argv[1] == "-v":
    print("1.0.0")
    sys.exit()
else:
    in_path = sys.argv[1]

in_path = sys.argv[1]
in_data = gpf.ll(in_path)

for line in in_data:
    split_line = line.split("\t")
    qseqid = split_line[0]
    sseqid = split_line[1]
    pident = split_line[2]
    length = split_line[3]
    mismatch = split_line[4]
    gapopen = split_line[5]
    qstart = split_line[6]
    qend = split_line[7]
    sstart = split_line[8]
    send = split_line[9]
    evalue = split_line[10]
    bitscore = split_line[11]
    staxids = split_line[12]
    sscinames = split_line[13]
    sskingdoms = split_line[14]
    sphylums = split_line[15]
    salltitles = split_line[16]

    out_line = "\t".join(
        [
            qseqid,
            staxids,
            bitscore,
            qseqid,
            sseqid,
            pident,
            length,
            mismatch,
            gapopen,
            qstart,
            qend,
            sstart,
            send,
            evalue,
            bitscore,
        ]
    )
    print(out_line)
