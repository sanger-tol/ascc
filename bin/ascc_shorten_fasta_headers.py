#!/usr/bin/env python3
VERSION = "1.1.0"
DESCRIPTION = f"""
---
Script for shortening FASTA headers, by splitting the header and keeping only the first element
Version: {VERSION}
---

Written by Eerik Aunin (ea10)

Modified by Damon-Lee Pointon (@dp24/@DLBPointon)

"""

# MIT License
#
# Copyright (c) 2020-2022 Genome Research Ltd.
#
# Author: Eerik Aunin (ea10@sanger.ac.uk)
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


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="ascc_shorten_fasta_headers",
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
    parser.add_argument("-v", "--version", action="version", version=VERSION)
    return parser.parse_args(argv)


def main(fasta_path, delimiter, allow_duplicate_headers):

    with tempfile.TemporaryDirectory() as tmp_dir:
        input_file = fasta_path
        if fasta_path.endswith(".gz") or fasta_path.endswith('.gz"'):
            input_file = "{}/input_file.fa".format(tmp_dir)
            gpf.run_system_command("gunzip -c {} > {}".format(fasta_path, input_file))

        headers_list = list()
        in_data = gpf.ll(input_file)

        for line in in_data:
            out_line = line
            if line.startswith(">"):
                if delimiter == "":
                    out_line = line.split()[0]
                else:
                    out_line = line.split(delimiter)[0]
                if out_line in headers_list and allow_duplicate_headers is False:
                    sys.stderr.write(
                        "Duplicate FASTA headers ({}) were found in the input file ({}) after truncating the headers with a delimiter\n".format(
                            out_line[1 : len(out_line)], fasta_path
                        )
                    )
                    sys.exit(1)
                headers_list.append(out_line)
            print(out_line)


if __name__ == "__main__":
    args = parse_args()
    main(args.fasta_path, args.delimiter, args.allow_duplicate_headers)
