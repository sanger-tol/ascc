#!/usr/bin/env python3

"""
Script to combine kmer count dimensionality reduction embedding csv files.
Adapted from Eerik Aunin (ea10@sanger.ac.uk), by William Eagles (we3@sanger.ac.uk)
"""

import pandas as pd
import argparse
from functools import reduce

def main(embedding_paths, output_name):

    # Load csv as dataframes.
    dataframe_list = list()
    for embedding_path in embedding_paths:
        df = pd.read_csv(f"{embedding_path}/kmers_dim_reduction_embeddings.csv")
        dataframe_list.append(df)

    # Merge DataFrames in list
    out_df = reduce(lambda left, right:     
                     pd.merge(left , right,
                              on = ["scaff"],
                              how = "outer"),
                     dataframe_list)

    # Output as csv.
    out_df.to_csv(output_name, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-v", "--version", action="version", version="1.0")
    parser.add_argument("-i", "--kmer_counts_file", nargs='+', type=str, help="List of paths to input CSV file with kmer counts")
    parser.add_argument("-o",  "--output_name", type=str, help="Path to folder where output file will be written")
    args = parser.parse_args()
    main(args.kmer_counts_file, args.output_name)