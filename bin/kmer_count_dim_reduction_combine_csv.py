#!/usr/bin/env python3

"""
Script to combine kmer count dimensionality reduction embedding csv files.
Adapted from Eerik Aunin (eeaunin@gmail.com), by William Eagles (we3@sanger.ac.uk)

This script takes CSV files containing dimensionality reduction embeddings and combines them
into a single CSV file. Each input CSV file should contain a 'scaff' column that will be used
for merging.
"""

import pandas as pd
import argparse
from functools import reduce
import sys
import os
import traceback


def main(embedding_paths, output_name):
    # Load csv as dataframes.
    dataframe_list = list()

    # Print debug information
    print(f"Processing {len(embedding_paths)} input paths", file=sys.stderr)
    for i, path in enumerate(embedding_paths):
        print(f"Input path {i+1}: {path} (type: {type(path).__name__})", file=sys.stderr)
        # Check if the path exists
        if os.path.exists(path):
            print(f"  Path exists", file=sys.stderr)
            if os.path.isdir(path):
                print(f"  Path is a directory", file=sys.stderr)
                # List contents of directory
                print(f"  Directory contents: {os.listdir(path)}", file=sys.stderr)
            elif os.path.isfile(path):
                print(f"  Path is a file", file=sys.stderr)
                # Get file size
                print(f"  File size: {os.path.getsize(path)} bytes", file=sys.stderr)
        else:
            print(f"  Path does not exist", file=sys.stderr)

    # First pass: just check all paths and print information
    for i, embedding_path in enumerate(embedding_paths):
        print(f"First pass - checking path {i+1}: {embedding_path}", file=sys.stderr)

        # Handle the case where embedding_path is a directory
        if os.path.isdir(embedding_path):
            # Look for the CSV file in the directory
            csv_path = os.path.join(embedding_path, "kmers_dim_reduction_embeddings.csv")
            if os.path.exists(csv_path):
                print(f"  Found CSV in directory: {csv_path}", file=sys.stderr)
                print(f"  File size: {os.path.getsize(csv_path)} bytes", file=sys.stderr)
            else:
                # Fallback to searching for any CSV file in the directory
                csv_files = [f for f in os.listdir(embedding_path) if f.endswith(".csv")]
                if csv_files:
                    found_csv = os.path.join(embedding_path, csv_files[0])
                    print(f"  Found CSV by searching directory: {found_csv}", file=sys.stderr)
                    print(f"  File size: {os.path.getsize(found_csv)} bytes", file=sys.stderr)
                else:
                    print(f"  Warning: No CSV files found in directory: {embedding_path}", file=sys.stderr)

        # Check if the path is a file
        elif os.path.isfile(embedding_path):
            print(f"  Path is a file", file=sys.stderr)
            print(f"  File size: {os.path.getsize(embedding_path)} bytes", file=sys.stderr)
        else:
            print(f"  Warning: {embedding_path} is neither a file nor a directory", file=sys.stderr)

    # Second pass: actually process the files
    for i, embedding_path in enumerate(embedding_paths):
        print(f"Second pass - processing path {i+1}: {embedding_path}", file=sys.stderr)

        # Handle the case where embedding_path is a directory
        if os.path.isdir(embedding_path):
            # Look for the CSV file in the directory
            csv_path = os.path.join(embedding_path, "kmers_dim_reduction_embeddings.csv")
            if os.path.exists(csv_path):
                embedding_path = csv_path
                print(f"  Using CSV in directory: {csv_path}", file=sys.stderr)
            else:
                # Fallback to searching for any CSV file in the directory
                csv_files = [f for f in os.listdir(embedding_path) if f.endswith(".csv")]
                if csv_files:
                    embedding_path = os.path.join(embedding_path, csv_files[0])
                    print(f"  Using CSV found by searching directory: {embedding_path}", file=sys.stderr)
                else:
                    print(f"  Warning: No CSV files found in directory: {embedding_path}", file=sys.stderr)
                    continue  # Skip this path and continue with the next one

        # Skip if the path doesn't exist or isn't a file
        if not os.path.isfile(embedding_path):
            print(f"  Warning: {embedding_path} is not a file, skipping", file=sys.stderr)
            continue

        try:
            print(f"  Reading CSV file: {embedding_path}", file=sys.stderr)
            df = pd.read_csv(embedding_path)
            print(f"  CSV columns: {df.columns.tolist()}", file=sys.stderr)

            # Verify the dataframe has a 'scaff' column
            if "scaff" not in df.columns:
                print(f"  Warning: {embedding_path} does not have a 'scaff' column, skipping", file=sys.stderr)
                continue

            print(
                f"  Successfully loaded {embedding_path} with {len(df)} rows and {len(df.columns)} columns",
                file=sys.stderr,
            )
            dataframe_list.append(df)
        except Exception as e:
            print(f"  Error reading {embedding_path}: {str(e)}", file=sys.stderr)
            # Print the traceback for more detailed error information
            import traceback

            traceback.print_exc(file=sys.stderr)
            continue

    # Check if we have any dataframes to merge
    if not dataframe_list:
        print("No valid CSV files were found to combine", file=sys.stderr)
        # Create an empty output file to satisfy Nextflow
        with open(output_name, "w") as f:
            f.write("scaff\n")  # Write header only
        print(f"Created empty output file: {output_name}", file=sys.stderr)
        return

    # Merge DataFrames in list
    out_df = reduce(lambda left, right: pd.merge(left, right, on=["scaff"], how="outer"), dataframe_list)

    # Output as csv.
    out_df.to_csv(output_name, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-v", "--version", action="version", version="1.0")
    parser.add_argument(
        "-i",
        "--kmer_counts_file",
        nargs="+",
        type=str,
        help="List of paths to input CSV files with dimensionality reduction embeddings",
    )
    parser.add_argument(
        "-o", "--output_name", type=str, help="Path to output CSV file where combined embeddings will be written"
    )
    args = parser.parse_args()
    main(args.kmer_counts_file, args.output_name)
