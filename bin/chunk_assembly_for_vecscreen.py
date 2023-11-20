#!/usr/bin/env python3
"""
Script for chunking an assembly before running NCBI VecScreen. Adapted from a script by James Torrance.
The script was further refactored by Eerik Aunin and Yumi Sims
"""

from Bio import SeqIO
import argparse
import os


def process_record(record, chunk_num, threshold_length, overlap_length):
    """
    ID chunks within a certain threshold, this sequence is then excised
    If it falls outside of the calculated threshold then rest of sequence is output
    """
    if chunk_num * threshold_length < len(record) - (threshold_length + overlap_length):
        return record[chunk_num * threshold_length : ((chunk_num + 1) * threshold_length + overlap_length)]
    else:
        return record[chunk_num * threshold_length :]


def generate_records_to_write(record, threshold_length, overlap_length):
    """
    For record in input process the record
    generate new id of id + counter
    """
    counter = 0
    for i in range((len(record) - 1) // threshold_length + 1):
        record_slice = process_record(record, i, threshold_length, overlap_length)
        counter += 1
        record_slice.id = record_slice.id + ".chunk_{}".format(counter)
        return record_slice


def main(fasta_input_file, fasta_output_file, threshold_length):
    fasta_input_file = os.path.abspath(fasta_input_file)
    fasta_output_file = os.path.abspath(fasta_output_file)

    overlap_length = int(threshold_length / 10)
    minimum_record_size = 11

    try:
        with open(fasta_output_file, "w") as fasta_output_handle:
            for record in SeqIO.parse(fasta_input_file, "fasta"):
                if len(record) >= minimum_record_size:
                    records_to_write = generate_records_to_write(record, threshold_length, overlap_length)
                    SeqIO.write(records_to_write, fasta_output_handle, "fasta")

    except Exception as e:
        print("An error occurred: {}".format(e))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta_input_file", type=str, help="Path to input FASTA file")
    parser.add_argument("fasta_output_file", type=str, help="Path for FASTA output file")
    parser.add_argument("--threshold", type=int, default=500000, help="Threshold length of sequence")
    parser.add_argument("-v", "--version", action="version", version="1.0")
    args = parser.parse_args()
    main(args.fasta_input_file, args.fasta_output_file, args.threshold)
