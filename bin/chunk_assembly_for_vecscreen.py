#!/usr/bin/env python3
"""
Script for chunking an assembly before running NCBI VecScreen. Adapted from a script by James Torrance.
The script was further refactored by Eerik Aunin and Yumi Sims
"""

from Bio import SeqIO
import argparse
import os


def process_record(record, chunk_num, threshold_length, overlap_length):
    return (
        record[chunk_num * threshold_length : ((chunk_num + 1) * threshold_length + overlap_length)]
        if chunk_num * threshold_length < len(record) - (threshold_length + overlap_length)
        else record[chunk_num * threshold_length :]
    )


def generate_records_to_write(record, threshold_length, overlap_length):
    return [
        (
            record_slice.id + ".chunk_{}".format(chunk_num + 1),
            "",
            record_slice,
        )
        for chunk_num, record_slice in enumerate(
            process_record(record, i, threshold_length, overlap_length)
            for i in range((len(record) - 1) // threshold_length + 1)
        )
    ]


def main(fasta_input_file, fasta_output_file):
    fasta_input_file = os.path.abspath(fasta_input_file)
    fasta_output_file = os.path.abspath(fasta_output_file)

    threshold_length = 500000
    overlap_length = int(threshold_length / 10)
    minimum_record_size = 11

    try:
        with open(fasta_output_file, "w") as fasta_output_handle:
            for record in SeqIO.parse(fasta_input_file, "fasta"):
                if len(record) >= minimum_record_size:
                    records_to_write = generate_records_to_write(record, threshold_length, overlap_length)
                    SeqIO.write(records_to_write, fasta_output_handle, "fasta")

    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta_input_file", type=str, help="Path to input FASTA file")
    parser.add_argument("fasta_output_file", type=str, help="Path for FASTA output file")
    parser.add_argument("-v", "--version", action="version", version="1.0")
    args = parser.parse_args()
    main(args.fasta_input_file, args.fasta_output_file)
