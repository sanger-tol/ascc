#!/usr/bin/env python3

import os
import argparse

"""
A simple script to generate csv file

Written by Damon-Lee Pointon (dp24/DLBPointon)
"""


def parse_args():
    parser = argparse.ArgumentParser(description="Generate a csv file for BTK")
    parser.add_argument("sample_name", type=str, help="Name of sample")
    parser.add_argument("pacbio_path", type=str, help="Path containing the pacbio files")
    parser.add_argument("-v", "--version", action="version", version="1.0.0")
    return parser.parse_args()


def main():
    args = parse_args()

    data_list = []

    data_list.append("sample,datatype,datafile\n")
    for file in os.listdir(args.pacbio_path):
        if file.endswith(".fasta.gz"):
            data_list.append(f"{args.sample_name},pacbio,{args.pacbio_path}/{file}\n")

    with open("samplesheet.csv", "w") as file:
        file.write("".join(data_list))


if __name__ == "__main__":
    main()
