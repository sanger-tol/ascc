#!/usr/bin/env python3
"""
Script for calculating coordinates in whole sequence from coordinates of sequence chunks in BLAST result files.
The BLAST file format is what is required by BlobToolKit

Written by Eerik Aunin @eeaunin

Adapted by Damon-Lee Pointon @DLBPointon
"""

import general_purpose_functions as gpf
import argparse
import sys


def process_nucleotide_blast_file(in_path):
    """
    Converts coordinates in a nucleotide BLAST file
    #these are the columns: qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    #the second column is taxid
    #blast coordinates are 1 based (http://alternateallele.blogspot.com/2012/03/genome-coordinate-cheat-sheet.html)
    """

    in_data = gpf.l(in_path)
    for line in in_data:
        try:
            split_line = line.split("\t")
            if len(split_line) < 11:  # Skip malformed lines
                sys.stderr.write(f"Warning: Skipping malformed line: {line}\n")
                continue

            field1_split = split_line[0].split("_sliding:")
            qseqid = field1_split[0]
            qoffset = int(field1_split[1].split("-")[0])
            qstart = int(split_line[9])
            qend = int(split_line[10])
            new_qstart = qstart + qoffset
            new_qend = qend + qoffset
            split_line[0] = qseqid
            split_line[3] = qseqid
            split_line[9] = str(new_qstart)
            split_line[10] = str(new_qend)
            out_line = "\t".join(split_line)
            print(out_line)
        except Exception as e:
            sys.stderr.write(f"Error processing line: {line}\nError: {str(e)}\n")
            continue


def process_diamond_file(in_path):
    """
    Converts coordinates in a Diamond BLASTX output file
    #outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms sphylums salltitles
    """

    in_data = gpf.l(in_path)
    for line in in_data:
        try:
            split_line = line.split("\t")
            if len(split_line) < 8:  # Skip malformed lines
                sys.stderr.write(f"Warning: Skipping malformed line: {line}\n")
                continue

            field1_split = split_line[0].split("_sliding:")
            qseqid = field1_split[0]
            qoffset = int(field1_split[1].split("-")[0])

            qstart = int(split_line[6])
            qend = int(split_line[7])

            new_qstart = qstart + qoffset
            new_qend = qend + qoffset

            split_line[0] = qseqid
            split_line[6] = str(new_qstart)
            split_line[7] = str(new_qend)
            out_line = "\t".join(split_line)
            print(out_line)
        except Exception as e:
            sys.stderr.write(f"Error processing line: {line}\nError: {str(e)}\n")
            continue


def main(in_path, blast_type):
    if blast_type == "nucleotide":
        process_nucleotide_blast_file(in_path)
    elif blast_type == "diamond":
        process_diamond_file(in_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-v", action="version", version="1.0")
    parser.add_argument("in_path", type=str, help="Path to BLAST results file")
    parser.add_argument(
        "blast_type", type=str, help="BLAST type: 'nucleotide' or 'diamond'", choices=["nucleotide", "diamond"]
    )
    args = parser.parse_args()
    main(args.in_path, args.blast_type)
