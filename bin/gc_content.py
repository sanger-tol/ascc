#!/usr/bin/env python3
"""
Script for finding the GC content of each sequence in a multiFASTA file

Originally Written by Eerik Aunin @eeaunin
Refactored by Yumi Sims @yumisims
Adapted by Damon-Lee Pointon @DLBPointon
"""

import argparse
import general_purpose_functions as gpf


def main(fasta_path):
    fasta_data = gpf.read_fasta_in_chunks(fasta_path)
    results = [
        (header.split()[0], "{:.6f}".format((seq.upper().count("G") + seq.upper().count("C")) / max(1, len(seq))))
        for header, seq in fasta_data
    ]
    for header, gc_content_string in results:
        print(f"{header}\t{gc_content_string}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta_path", type=str, help="Path to input FASTA file")
    parser.add_argument("-v", action="version", version="1.0")
    args = parser.parse_args()
    main(args.fasta_path)
