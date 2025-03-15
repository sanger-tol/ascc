#!/usr/bin/env python3
"""
Script for converting BLAST results to a format suitable for BlobToolKit
Input: BLAST results file in BTK format (qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue)
Output: BLAST hits file for BlobToolKit (qseqid staxids bitscore)

This script takes a BLAST results file in BTK format and extracts just the first three columns
that are required by BlobToolKit for visualization.
"""

import sys
import argparse

def convert_blast_to_hits(input_file, output_file):
    """
    Convert BLAST results to BlobToolKit hits format
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                # Extract just the first three columns: qseqid, staxids, bitscore
                outfile.write(f"{fields[0]}\t{fields[1]}\t{fields[2]}\n")
            else:
                sys.stderr.write(f"Warning: Skipping line with unexpected format: {line}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert BLAST results to BlobToolKit hits format")
    parser.add_argument("--input", required=True, help="Input BLAST results file in BTK format")
    parser.add_argument("--output", required=True, help="Output file in BlobToolKit hits format")
    
    args = parser.parse_args()
    
    convert_blast_to_hits(args.input, args.output)
    
    sys.stderr.write(f"Converted BLAST hits written to {args.output}\n")
