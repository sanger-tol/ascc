#!/usr/bin/env python3
"""
Script for chunking an assembly before running NCBI VecScreen. Adapted from a script by James Torrance, edited by Eerik Aunin and Yumi Sims
"""

from Bio import SeqIO
import argparse
import os


def process_record(record, threshold_length, overlap_length):
    """
    ID chunks within a certain threshold, this sequence is then excised
    If it falls outside of the calculated threshold then rest of sequence is output
    if a scaffold is shorter than or equal to the specified threshold_length, it will not be chunked, and the chunk ID won't be attached.
    """

    record_slices = [record[i : i + threshold_length + overlap_length] for i in range(0, len(record), threshold_length)]

    chunks = []
    for slice_count, record_slice in enumerate(record_slices, start=1):
        if len(record) > threshold_length:
            record_slice.id = "{}.chunk_{}".format(record_slice.id, slice_count)
        record_slice.description = ""
        chunks.append(record_slice)

    return chunks


def generate_records_to_write(record, threshold_length, overlap_length):
    """
    For record in input process the record
    generate new id of id + counter
    """
    return [
        (
            "{}.chunk_{}".format(record.id, chunk_num + 1),
            "",
            record_slice,
        )
        for chunk_num, record_slice in enumerate(process_record(record, threshold_length, overlap_length), start=1)
    ]


def main(fasta_input_file, fasta_output_file, threshold_length):
    fasta_input_file = os.path.abspath(fasta_input_file)
    fasta_output_file = os.path.abspath(fasta_output_file)

    overlap_length = int(threshold_length / 10)
    minimum_record_size = 11

    try:
        with open(fasta_input_file, "r") as fasta_input_handle, open(fasta_output_file, "w") as fasta_output_handle:
            for record in SeqIO.parse(fasta_input_handle, "fasta"):
                if len(record) >= minimum_record_size:
                    records_to_write = generate_records_to_write(record, threshold_length, overlap_length)
                    SeqIO.write([chunk[2] for chunk in records_to_write], fasta_output_handle, "fasta")

    except Exception as e:
        print("An error occurred: {}".format(e))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta_input_file", type=str, help="Path to input FASTA file")
    parser.add_argument("fasta_output_file", type=str, help="Path for FASTA output file")
    parser.add_argument("--threshold_length", type=int, default=500000, help="Threshold length for chunking")
    parser.add_argument("-v", "--version", action="version", version="1.0")
    args = parser.parse_args()
    main(args.fasta_input_file, args.fasta_output_file, args.threshold_length)
