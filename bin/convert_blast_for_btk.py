#!/usr/bin/env python3
"""
Script for converting BLAST results to a format suitable for BlobToolKit
Input: BLAST results file (format 6)
Output: BLAST results file in BlobToolKit format

This script takes a BLAST results file and reformats it to be compatible with BlobToolKit.
The input format is: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
The output format is: qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue
"""

import sys
import argparse

def convert_blast_for_btk(input_file, output_file):
    """
    Convert BLAST results to BlobToolKit format
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            if len(fields) >= 12:
                # Extract the required fields
                qseqid = fields[0]
                sseqid = fields[1]
                pident = fields[2]
                length = fields[3]
                mismatch = fields[4]
                gapopen = fields[5]
                qstart = fields[6]
                qend = fields[7]
                sstart = fields[8]
                send = fields[9]
                evalue = fields[10]
                bitscore = fields[11]
                
                # Create a new line with the format expected by BlobToolKit
                # qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue
                new_line = f"{qseqid}\t0\t{bitscore}\t{qseqid}\t{sseqid}\t{pident}\t{length}\t{mismatch}\t{gapopen}\t{qstart}\t{qend}\t{sstart}\t{send}\t{evalue}"
                outfile.write(new_line + '\n')
            else:
                sys.stderr.write(f"Warning: Skipping line with unexpected format: {line}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert BLAST results to BlobToolKit format")
    parser.add_argument("--input", required=True, help="Input BLAST results file")
    parser.add_argument("--output", required=True, help="Output file in BlobToolKit format")
    
    args = parser.parse_args()
    
    convert_blast_for_btk(args.input, args.output)
    
    sys.stderr.write(f"Converted BLAST results written to {args.output}\n")
