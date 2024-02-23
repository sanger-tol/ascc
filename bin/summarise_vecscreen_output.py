#!/usr/bin/env python3
"""
Script for reformatting output from VecScreen that has been processed using the VSlistTo1HitPerLine.awk script. 
This script removes no-hit entries from the report file and converts coordinates in assembly chunks to coordinates in the full assembly.
The original script was written in Perl by James Torrance. Rewritten by Eerik Aunin.
"""

import argparse
import os
import sys


def main(vecscreen_file, chunk_size):
    if os.path.isfile(vecscreen_file) is False:
        sys.stderr.write(f"The input file for this script ({vecscreen_file}) was not found\n")
        sys.exit(1)

    with open(vecscreen_file) as f:
        vecscreen_data = f.readlines()

    output_lines = list()

    for line in vecscreen_data:
        if line.startswith("VecScreen_No_Hits") is True:
            continue
        split_line = line.split()
        assert len(split_line) == 4
        vecscreen_result, seq_name, seq_start, seq_end = split_line
        seq_start = int(seq_start)
        seq_end = int(seq_end)
        seq_base_name = seq_name
        if ".chunk_" in seq_name:
            split_seq_name = seq_name.split(".chunk_")
            assert len(split_seq_name) == 2
            seq_base_name = split_seq_name[0]
            chunk_suffix = int(split_seq_name[1])
            seq_start += chunk_size * (chunk_suffix - 1)
            seq_end += chunk_size * (chunk_suffix - 1)
        corrected_line = "\t".join((vecscreen_result, seq_base_name, str(seq_start), str(seq_end)))
        if corrected_line not in output_lines:
            output_lines.append(corrected_line)

    for line in output_lines:
        print(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "vecscreen_file",
        type=str,
        help="Path to output file of VecScreen (run with -f3 flag), filtered with VSlistTo1HitPerLine.awk",
    )
    parser.add_argument(
        "--chunk_size",
        type=int,
        help="Chunk size of the chunked FASTA file that VecScreen was run with, in bp. Default: 500000",
        default=50000,
    )
    parser.add_argument("-v", "--version", action="version", version="1.0")
    args = parser.parse_args()
    main(args.vecscreen_file, args.chunk_size)
