#!/usr/bin/env python3

VERSION = "1.1.1"
DESCRIPTION = f"""
---
Script for checking if the TaxID given by the user exists in the NCBI taxdump
Version: {VERSION}
---

Written by Eerik Aunin (ea10)

Modified by Damon-Lee Pointon (@dp24/@DLBPointon)
Modified to add F-Strings 24/09/26

"""

import argparse
import os
import sys
import textwrap

import general_purpose_functions as gpf


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="find_taxid_in_taxdump",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(DESCRIPTION),
    )
    parser.add_argument("query_taxid", type=int, help="Query taxonomy ID")
    parser.add_argument(
        "taxdump_nodes_path",
        type=str,
        help="Path to the nodes.dmp file of NCBI taxdump (downloaded from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)",
    )
    parser.add_argument("-v", "--version", action="version", version=VERSION)

    return parser.parse_args(argv)


def main(query_taxid, taxdump_nodes_path):
    if query_taxid == -1:
        sys.exit(0)
    query_taxid = str(query_taxid)
    if os.path.isfile(taxdump_nodes_path) is False:
        sys.stderr.write(f"The NCBI taxdump nodes file ({taxdump_nodes_path}) was not found\n")
        sys.exit(1)
    nodes_data = gpf.ll(taxdump_nodes_path)
    taxid_found_flag = False
    for counter, line in enumerate(nodes_data):
        split_line = line.split("|")
        if len(split_line) > 2:
            taxid = split_line[0].strip()
            if taxid == query_taxid:
                taxid_found_flag = True
                break
        else:
            sys.stderr.write(
                f"Failed to parse the NCBI taxdump nodes.dmp file ({taxdump_nodes_path}) at line {counter + 1}:\n"
            )
            sys.stderr.write(line + "\n")
            sys.exit(1)

    if taxid_found_flag is False:
        sys.stderr.write(
            f"The TaxID given by the user ({query_taxid}) was not found in the NCBI taxdump nodes.dmp file ({taxdump_nodes_path})\n"
        )
        sys.exit(1)


if __name__ == "__main__":
    args = parse_args()
    main(args.query_taxid, args.taxdump_nodes_path)
