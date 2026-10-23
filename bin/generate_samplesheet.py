#!/usr/bin/env python3

import os
import sys
import argparse

"""
A simple script to generate a csv file required for the sanger-tol/blobtoolkit pipeline-module.

Required input include the sample ID and the mapped BAM file generated with PacBio data and input FASTA assembly

Written by Damon-Lee Pointon (dp24/DLBPointon)
"""


def parse_args():
    parser = argparse.ArgumentParser(description="Generate a csv file for BTK")
    parser.add_argument("sample_name", type=str, help="Name of sample")
    parser.add_argument(
        "path_to_reads",
        type=str,
        help="Path containing the PacBio reads",
    )
    parser.add_argument(
        "reads_layout", type=str, help="Whether the reads are SINGLE or PAIRED end"
    )
    parser.add_argument("-v", "--version", action="version", version="1.2.0")

    return parser.parse_args()


def main():
    args = parse_args()

    data_list = []

    data_list.append("sample,datatype,datafile,library_layout\n")

    [
        data_list.append(
            f"{args.sample_name},pacbio,{args.path_to_reads}{file},{args.reads_layout}\n"
        )
        for file in os.listdir(args.path_to_reads)
        if file.endswith(".fasta.gz") or file.endswith(".fa.gz")
    ]

    if len(data_list) <= 1:
        sys.exit("I was expecting at least one FASTA.GZ file")

    with open("samplesheet.csv", "w") as file:
        file.write("".join(data_list))


if __name__ == "__main__":
    main()
