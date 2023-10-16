#!/usr/bin/env python3

"""
Pacbio Barcode Check
------------------------
Looks for Pacbio barcodes in ref and data.

If any User-supplied barcodes arn't in data then pipeline dies.
If no barcodes are found then pipeline carries on

Based on a standard operating procedure developed by James Torrance

Originally written by Eerik Aunin @eeaunin

Adapted by Damon-Lee Pointon @DLBPointon
"""

from pathlib import Path
import os
import sys
import argparse
import general_purpose_functions as gpf


def detect_barcodes_from_read_file_names(barcodes_fasta_path, pacbio_read_files):
    """
    Reads barcode names, e.g. bc1001_BAK8A_OA from the barcodes FASTA file headers and then looks for these names in PacBio read file names
    Returns a list of barcode names detected in PacBio reads file headers
    """
    barcodes_fasta_data = gpf.l(barcodes_fasta_path)
    barcode_names = [n.split(">")[1] for n in barcodes_fasta_data if n.startswith(">")]
    if len(barcode_names) == 0:
        sys.stderr.write(
            f"Failed to read PacBio multiplexing barcode names from the specified file {barcodes_fasta_data}\n"
        )
        sys.exit(1)
    detected_barcodes = list()
    for barcode_name in barcode_names:
        for read_file in pacbio_read_files:
            query_string = f".ccs.{barcode_name}--{barcode_name}."
            if query_string in read_file:
                detected_barcodes.append(barcode_name)
    return detected_barcodes


def check_if_barcodes_exist_in_barcodes_fasta(barcodes_list, barcodes_fasta_path):
    """
    Checks if the specified barcodes exist in the barcode sequences FASTA file, exits with an error message if a barcode is not found
    """
    barcodes_fasta_data = gpf.l(barcodes_fasta_path)
    barcode_names_in_fasta = [n.split(">")[1] for n in barcodes_fasta_data if n.startswith(">")]
    for barcode in barcodes_list:
        if barcode not in barcode_names_in_fasta:
            sys.stderr.write(
                f"The PacBio multiplexing barcode ({barcode}) was not found in the barcode sequences file ({barcodes_fasta_path})\n"
            )
            sys.exit(1)

    # If this print statement is reached, all user-supplied codes are present.
    print("The query barcodes exist in the barcodes database file")


def main(barcodes_fasta_path, pacbio_read_files, pacbio_multiplexing_barcode_names):
    pacbio_read_files = pacbio_read_files.split(",")

    barcodes_list = []
    if len(pacbio_multiplexing_barcode_names) > 0:
        barcodes_list = pacbio_multiplexing_barcode_names.strip("[").strip("]").split(",")

    if os.path.isfile(barcodes_fasta_path) is False:
        sys.stderr.write(
            "FASTA file with PacBio multiplexing barcode sequences ({barcodes_fasta_path}) was not found\n"
        )
        sys.exit(1)

    if len(barcodes_list) == 0:
        barcodes_list = detect_barcodes_from_read_file_names(barcodes_fasta_path, pacbio_read_files)

    # Here script should break successfully
    if len(barcodes_list) == 0:
        sys.stderr.write(
            "Skipping the PacBio barcodes check, as no barcodes were specified by the user and no barcodes were found in PacBio read file names\n"
        )
        print("Skipping the PacBio barcodes check, as no barcodes were specified by the user and no barcodes were found in PacBio read file names\n")
        sys.exit(0)

    check_if_barcodes_exist_in_barcodes_fasta(barcodes_list, barcodes_fasta_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-b", "--barcode_fasta", type=str, help="Pacbio Barcode FASTA file")

    parser.add_argument("-p", "--pacbio_reads", type=str, help="Pacbio Read FASTA.gz files")

    parser.add_argument("-m", "--multiplex_name", type=str, help="Pacbio Multiplex Barcode Name")

    parser.add_argument("-v", action="version", version="1.0")
    args = parser.parse_args()
    main(args.barcode_fasta, args.pacbio_reads, args.multiplex_name)

    # DOWNSTREAM
    #   barcodes_fasta_path  passed into this script is also required for make_blast_db
    #   stdout(barcodes_list).split_csv().set(i) -> filter_barcode(i)
