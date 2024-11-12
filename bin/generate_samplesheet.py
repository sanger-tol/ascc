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
        "mapped_bam_file",
        type=str,
        help="Path containing the mapped BAM generated with PacBio data and the ASCC input assembly",
    )
    parser.add_argument("-v", "--version", action="version", version="1.0.0")
    return parser.parse_args()


def main():
    args = parse_args()

    data_list = []

    data_list.append("sample,datatype,datafile\n")
    if args.mapped_bam_file.endswith(".bam"):
        data_list.append(f"{args.sample_name},pacbio,{args.mapped_bam_file}\n")
    else:
        sys.exit("I was expecting a mapped BAM file")

    with open(f"{args.sample_name}_samplesheet.csv", "w") as file:
        file.write("".join(data_list))


if __name__ == "__main__":
    main()
