#!/usr/bin/env python3
"""
Script for counting kmer frequencies per sequence in a FASTA file
Output (STDOUT): kmer counts as a CSV table
Developed by Eerik Aunin (ea10@sanger.ac.uk)
"""

import argparse
import general_purpose_functions as gpf
import kcounter
from collections import OrderedDict
import pandas as pd


def main(fasta_path, out_path, kmer_size):
    fasta_data = gpf.read_fasta_in_chunks(fasta_path)
    nucleotides_collection = list()
    for header, seq in fasta_data:
        seq = seq.upper()
        seq_len = len(seq)
        nucleotides_dict = kcounter.count_kmers(seq, kmer_size, canonical_kmers=True)
        relative_counts_dict = OrderedDict()
        relative_counts_dict["header"] = header
        relative_counts_dict["seq_len"] = seq_len
        for kmer in nucleotides_dict:
            kmer_relative_count = nucleotides_dict[kmer] / seq_len
            relative_counts_dict[kmer] = kmer_relative_count
        nucleotides_collection.append(relative_counts_dict)
    df = pd.DataFrame(nucleotides_collection)
    df = df.fillna(0)
    df.to_csv(out_path, index=False)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-v", "--version", action="version", version="1.0")
    parser.add_argument("fasta_path", type=str, help="Path to input FASTA file")
    parser.add_argument("out_path", type=str, help="Path for output CSV file")
    parser.add_argument("--kmer_size", type=int, help="kmer size (bp). Default: 7", default=7)
    args = parser.parse_args()
    main(args.fasta_path, args.out_path, args.kmer_size)
        