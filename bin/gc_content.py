#!/usr/bin/env python3
"""
Script for finding the GC content of each sequence in a multiFASTA file
"""

import argparse
import general_purpose_functions as gpf


def main(fasta_path):
    fasta_data = gpf.read_fasta_in_chunks(fasta_path)
    for header, seq in fasta_data:
        header = header.split()[0]
        seq = seq.upper()
        gc_content = None
        gc_count = seq.count("G") + seq.count("C")
        seq_len = len(seq)
        if seq_len > 0:
            gc_content = gc_count / seq_len
        gc_content_string = "{:.6f}".format(gc_content)
        print("{}\t{}".format(header, gc_content_string))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta_path", type=str, help="Path to input FASTA file")
    parser.add_argument("-v", action="version", version="1.0")
    args = parser.parse_args()
    main(args.fasta_path)
