#!/usr/bin/env python3
VERSION = "1.2.0"
DESCRIPTION = f"""
---
Script for sanitising FASTA headers and sequences:
- Shortens headers by splitting by whitespace and keeping only the first element
- Replaces problematic characters in headers (commas, spaces, etc.) with underscores
- Converts sequences to uppercase and replaces non-ATGC bases with N
Version: {VERSION}
---

Written by Eerik Aunin (ea10)
Modified by Damon-Lee Pointon (@dp24/@DLBPointon)
Further modified by Eerik Aunin (ea10)

"""

# MIT License
#
# Copyright (c) 2020-2022 Genome Research Ltd.
#
# Author: Eerik Aunin (eeaunin@gmail.com)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import general_purpose_functions as gpf
import argparse
import textwrap
import sys
import tempfile
import re


def is_all_n_sequence(seq):
    """Return True if sequence consists entirely of N's."""
    return all(base == "N" for base in seq.strip().upper())


def sanitise_sequence(seq):
    """Convert sequence to uppercase and replace any non-ATGC bases with N."""
    seq = seq.upper()
    return re.sub(r"[^ATGC]", "N", seq)


def sanitise_header(header):
    """Replace problematic characters in FASTA headers with underscores."""
    # Remove the '>' character if present at the start
    if header.startswith(">"):
        header = header[1:]

    # Replace problematic characters with underscores
    sanitised = re.sub(r"[,;\s|:]", "_", header)

    # Add back the '>' character
    return ">" + sanitised


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="sanitise_input_fasta_file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(DESCRIPTION),
    )
    parser.add_argument("fasta_path", type=str, help="Path to input FASTA file")
    parser.add_argument(
        "--delimiter",
        type=str,
        help="Delimiter string for splitting FASTA headers. Default: any whitespace character",
        default="",
    )
    parser.add_argument("--allow_duplicate_headers", dest="allow_duplicate_headers", action="store_true")
    parser.add_argument(
        "--keep_n_sequences", action="store_true", help="Keep sequences that are all Ns (default: False)"
    )
    parser.add_argument("-v", "--version", action="version", version=VERSION)
    return parser.parse_args(argv)


def main(fasta_path, delimiter, allow_duplicate_headers, keep_n_sequences=False):
    with tempfile.TemporaryDirectory() as tmp_dir:
        input_file = fasta_path
        if fasta_path.endswith(".gz") or fasta_path.endswith('.gz"'):
            input_file = "{}/input_file.fa".format(tmp_dir)
            gpf.run_system_command("gunzip -c {} > {}".format(fasta_path, input_file))

        headers_list = list()
        headers_with_commas = 0
        in_data = gpf.ll(input_file)

        current_header = None
        current_sequence = []

        def process_sequence():
            if current_header and current_sequence:
                sequence = "".join(current_sequence)
                if keep_n_sequences or not is_all_n_sequence(sequence):
                    print(current_header)
                    print(sequence)
                else:
                    sys.stderr.write("Skipping all-N sequence: {}\n".format(current_header[1:].strip()))

        for line in in_data:
            if line.startswith(">"):
                # Process previous sequence if it exists
                process_sequence()

                # Start new sequence
                original_header = line.strip()
                if delimiter == "":
                    current_header = original_header.split()[0]
                else:
                    current_header = original_header.split(delimiter)[0]

                # Check for commas in the original header
                if "," in original_header:
                    headers_with_commas += 1

                # Sanitise the header
                current_header = sanitise_header(current_header)

                if current_header in headers_list and allow_duplicate_headers is False:
                    sys.stderr.write(
                        "Duplicate FASTA headers ({}) were found in the input file ({}) after truncating the headers with a delimiter\n".format(
                            current_header[1:], fasta_path
                        )
                    )
                    sys.exit(1)
                headers_list.append(current_header)
                current_sequence = []
            else:
                # Add sanitised sequence line
                current_sequence.append(sanitise_sequence(line))

        # Process the last sequence
        process_sequence()

        # Print warning about headers with commas
        if headers_with_commas > 0:
            sys.stderr.write(
                "Warning: {} FASTA header(s) contained commas that were replaced with underscores\n".format(
                    headers_with_commas
                )
            )


if __name__ == "__main__":
    args = parse_args()
    main(args.fasta_path, args.delimiter, args.allow_duplicate_headers, args.keep_n_sequences)
