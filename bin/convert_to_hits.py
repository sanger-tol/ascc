#!/usr/bin/env python3

import general_purpose_functions as gpf
from collections import OrderedDict
import textwrap
import argparse

VERSION = "2.0.0"

DESCRIPTION = """
Script for getting the top hits of Diamond BLASTX against the nr database. Top hits per each scaffold are determined from the BLASTX results of scaffold chunks.
The output is reformatted into a CSV table
Argument1: path to Diamond BLASTX results, with the following columns:
"qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "end", "evalue", "bitscore", "staxids", "sscinames", "sskingdoms", "sphylums", "salltitles"

Version: {VERSION}
---
Written by Eerik Aunin

Re-Written by Damon-Lee Pointon (dp24/DLBPointon)
"""


def save_file(output_list, name, prefix):
    with open(f"{prefix}_{name}_diamond_blastx_top_hits.csv", "w") as f:
        for line in output_list:
            f.write(f"{line}\n")


def main(in_path, diamond_database_title, out_prefix):
    in_data = gpf.ll(in_path)

    top_hits_dict = OrderedDict()
    output = []

    colnames = (
        "scaff",
        "diamond_{0}_bitscore",
        "diamond_{0}_staxids",
        "diamond_{0}_sscinames",
        "diamond_{0}_sskingdoms",
        "diamond_{0}_sphylums",
        "diamond_{0}_salltitles",
    )
    colnames = [n.format(diamond_database_title) for n in colnames]

    for line in in_data:
        line = line.replace(",", " ")
        line = line.replace("<>", " ")
        split_line = line.split("\t")
        seq_name = split_line[0]
        seq_name = seq_name.split("_sliding")[0]
        bitscore = float(split_line[11])
        if seq_name not in top_hits_dict:
            top_hits_dict[seq_name] = dict()
            top_hits_dict[seq_name]["bitscore"] = bitscore
            top_hits_dict[seq_name]["line"] = line
        else:
            if bitscore > top_hits_dict[seq_name]["bitscore"]:
                top_hits_dict[seq_name]["bitscore"] = bitscore
                top_hits_dict[seq_name]["line"] = line

    output.append(",".join(colnames))
    for seq_name in top_hits_dict:
        out_line = top_hits_dict[seq_name]["line"]
        out_line = out_line.split("\t")
        out_line = ",".join(out_line[11 : len(out_line)])
        out_line = seq_name + "," + out_line
        output.append(out_line)

    save_file(output, diamond_database_title, out_prefix)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="Convert File to hits",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(DESCRIPTION),
    )
    parser.add_argument("in_path", type=str, help="Path to Diamond BLASTX results")
    parser.add_argument(
        "diamond_database_title",
        type=str,
        help="Name of the Diamond database (e.g. nr or Uniprot)",
    )
    parser.add_argument("-v", "--version", action="version", version=VERSION)
    parser.add_argument("-o", "--out_prefix", type=str, help="Output Prefix")
    args = parser.parse_args()
    main(args.in_path, args.diamond_database_title, args.out_prefix)
