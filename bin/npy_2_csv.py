#!/usr/bin/env python3

import sys
import textwrap
import argparse
import numpy as np
from itertools import product

VERSION = "1.0.0"
DESCRIPTION = f"""
NPY-2-CSV
Version = {VERSION}
---------
Written by Damon-Lee Pointon (DLBPointon, dp24), (14/05/2025)

NOTE: This script aims to be temporary until downstream scripts are
modified to take an input npy file.

This script converts an input NPY file from kmer-counter (Claudia Webber),
and converts it into a CSV file formatted for use in the downstream kmer
scripts written by Eerik Anuin.

This requires the KMER sequence, which is re-generated in the function `kmer_generator`.

"""


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="NPY to CSV",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(DESCRIPTION),
    )

    parser.add_argument(
        "-f",
        "--fasta",
        type=str,
        required=True,
        help="Path to the original fasta file used in kmer_counter.",
    )

    parser.add_argument(
        "-n",
        "--npy",
        type=str,
        required=True,
        help="Path to the output NPY of kmer-counter.",
    )

    parser.add_argument(
        "-k",
        "--kmer",
        type=int,
        default=7,
        help="The kmer length used in kmer-counter.",
    )

    parser.add_argument(
        "-o", "--output", type=str, default="npy.csv", help="The output CSV file naming"
    )

    parser.add_argument("-v", "--version", action="version", version=VERSION)

    return parser.parse_args(argv)


def kmer_generator(kmer_length):
    """
    Regenerate the kmer list used by the kmer_counter script
    based on function provided by the kmer_counter docs.
    """
    return ["".join(i) for i in product("ACGT", repeat=kmer_length)]


def return_canonical_only(kmer_list):
    """
    Return the first half of the kmer list this will contain only
    kmers not found in its reverse complement
    """
    return kmer_list[: (int(len(kmer_list) / 2))]


def get_complement(kmer) -> str:
    """
    Generate the complement of a given sequence
    """
    complement_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join([complement_dict[i] for i in kmer])


def complement_check(kmer_list) -> None:
    """
    Ensure that there are no complementary kmers found in the
    list.
    """
    for i in kmer_list:
        comp_i = get_complement(i)
        if comp_i in kmer_list:
            sys.exit("Non-canonical kmer found!!!")


def get_fasta_data(fasta) -> dict:
    """
    Open the fasta file and generate a dictionary of:
        {sequence_header: length_of_sequence}
    """
    data = {}
    tmp = ""
    line_counter = 0
    with open(fasta) as input:
        for line in input:
            if line.startswith(">") and line_counter == 0:
                tmp = line[1:].strip()
            elif line.startswith(">") and line_counter != 0:
                data[line[1:].strip()] = line_counter
                line_counter = 0
            else:
                line_counter += len(line.strip())
        data[tmp] = line_counter
    return data


def convert_dict_to_list(fasta_dict) -> list:
    """
    Convert the dictionary of sequence/length into a list of lists
    """
    return [[x, y] for x, y in fasta_dict.items()]


def main() -> None:
    args = parse_args()

    # kmer generator
    kmers = kmer_generator(args.kmer)
    canonical_kmers = return_canonical_only(kmers)

    # Validate kmer list
    complement_check(canonical_kmers)

    # New output header line
    header_line = ["header", "seq_len"] + canonical_kmers

    # Get fasta data
    fasta_dict = get_fasta_data(args.fasta)
    fasta_list_list = convert_dict_to_list(fasta_dict)

    # Read the npy file
    read = np.load(args.npy)

    # Set assertions for length and width of data
    assert len(read) == len(fasta_list_list)
    assert len(read[0]) == len(canonical_kmers)

    # Zip the fasta data with the npy data
    # both of which are ordered by order of the input fasta
    mix = zip(fasta_list_list, read.tolist())

    # open output file and start dumping data in the right format
    with open(args.output, "w") as output:
        output.write(",".join(header_line) + "\n")
        for line in mix:
            output_line = list(map(str, line[0])) + list(map(str, line[1]))

            # assert the length of (scaffold + counts) lines == the length of the header line
            assert len(output_line) == len(canonical_kmers) + 2
            output.write(",".join(output_line) + "\n")


if __name__ == "__main__":
    main()
