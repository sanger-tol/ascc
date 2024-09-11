#!/usr/bin/env python3

VERSION = "1.1.0"
DESCRIPTION = f"""
---
Script for filtering a FASTA file by sequence length. By default, sequences shorter than a cutoff value will be removed.
Version = {VERSION}
---

Written by Eerik Aunin (ea10)

Modified by Damon-Lee Pointon (@dp24/@DLBPointon)

"""

import general_purpose_functions as gpf
import textwrap
import argparse
import sys
import os


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="filter_fasta_by_length",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(DESCRIPTION),
    )
    parser.add_argument("in_path", type=str, help="Path to input FASTA file")
    parser.add_argument("cutoff", type=int, help="Cutoff value for filtering")
    parser.add_argument(
        "-l",
        "--low_pass",
        dest="low_pass",
        action="store_true",
        help="Optional: low pass filtering mode (sequences longer than the cutoff value will be removed)",
    )
    parser.add_argument(
        "--remove_original_fasta",
        action="store_true",
        help="Optional: remove the input FASTA file after creating the filtered FASTA file",
    )
    parser.add_argument("-v", "--version", action="version", version=VERSION)

    return parser.parse_args(argv)


def main(args):
    fasta_path = os.path.abspath(args.in_path)
    if (
        args.cutoff == -1
    ):  # When this script is used as a part of a pipeline, -1 can be assigned as a value for the cutoff to indicate that no filtering should be done
        sys.stderr.write(f"The input FASTA sequences ({fasta_path}) will not be filtered by length\n")
        # sys.exit(0)
    retained_seq_count = 0
    fasta_data = gpf.read_fasta_in_chunks(fasta_path)
    for header, seq in fasta_data:
        if args.cutoff == -1:
            print(">" + header)
            gpf.print_with_fixed_row_length(seq, 80)
        else:
            seq_len = len(seq)
            seq_len_ok_flag = True
            if args.low_pass == True:
                if seq_len > args.cutoff:
                    seq_len_ok_flag = False
                    sys.stderr.write(
                        f"Low pass filtering of FASTA sequences by length: removing sequence {header} from the assembly because its length ({seq_len}) exceeds the length cutoff ({args.cutoff})\n"
                    )
            else:
                if seq_len < args.cutoff:
                    seq_len_ok_flag = False
                    sys.stderr.write(
                        f"High pass filtering of FASTA sequences by length: removing sequence {header} from the assembly because its length ({seq_len}) is below the length cutoff ({args.cutoff})\n"
                    )
            if seq_len_ok_flag == True:
                retained_seq_count += 1
                print(">" + header)
                gpf.print_with_fixed_row_length(seq, 80)
                # print(header, seq_len, seq_len_ok_flag)
    if args.cutoff != -1 and retained_seq_count == 0:
        sys.stderr.write(
            f"No sequences remain in the FASTA file {fasta_path} after filtering the sequences by length (cutoff: {args.cutoff} bp)\n"
        )
        sys.exit(1)

    if args.remove_original_fasta is True:
        os.remove(fasta_path)


if __name__ == "__main__":
    main(parse_args())
