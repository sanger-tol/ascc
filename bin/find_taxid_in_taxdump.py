#!/usr/bin/env python3

VERSION = "1.1.0"
DESCRIPTION = f"""
---
Script for checking if the TaxID given by the user exists in the NCBI taxdump
Version: {VERSION}
---

Written by Eerik Aunin (ea10)

Modified by Damon-Lee Pointon (@dp24/@DLBPointon)

"""

import argparse
import textwrap
import general_purpose_functions as gpf
import sys
import os


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
        sys.stderr.write(
            "The NCBI taxdump nodes file ({}) was not found\n".format(
                taxdump_nodes_path
            )
        )
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
                "Failed to parse the NCBI taxdump nodes.dmp file ({}) at line {}:\n".format(
                    taxdump_nodes_path, counter + 1
                )
            )
            sys.stderr.write(line + "\n")
            sys.exit(1)

    if taxid_found_flag is False:
        sys.stderr.write(
            "The TaxID given by the user ({}) was not found in the NCBI taxdump nodes.dmp file ({})\n".format(
                query_taxid, taxdump_nodes_path
            )
        )
        sys.exit(1)


if __name__ == "__main__":
    args = parse_args()
    main(args.query_taxid, args.taxdump_nodes_path)
