#!/usr/bin/env python3
"""
Script for getting lineage for Kraken results
Developed by Eerik Aunin (eeaunin@gmail.com)
"""

import general_purpose_functions as gpf
import sys
from collections import OrderedDict
import pandas as pd
import os
import signal
import argparse


def load_kraken_results(kraken_results_path):
    """
    Reads sequence names and taxid values from Kraken results file into a dictionary
    """
    kraken_dict = OrderedDict()
    kraken_data = gpf.ll(kraken_results_path)
    for line in kraken_data:
        if "(taxid " in line:
            split_line = line.split()
            if len(split_line) >= 5:
                seq_name = split_line[1]
                taxid = gpf.spl(line, "(taxid ", ")")
                if seq_name in kraken_dict:
                    sys.stderr.write("Duplicate read names found in input ({})\n".format(seq_name))
                    os.kill(os.getpid(), signal.SIGINT)
                else:
                    kraken_dict[seq_name] = taxid
            else:
                sys.stderr.write("Failed to parse Kraken output file line:\n{}\n".format(line))
        else:
            sys.stderr.write("No taxid found in input file line:\n{}\n".format(line))
    return kraken_dict


def load_lineage(lineage_dump_path, kraken_db_name):
    """
    Reads lineage information from NCBI rankedlineage.dmp file (downloaded from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz) as a dictionary of dictionaries.
    Keys: taxid values. Values: a dictionary where the keys are taxonomic unit names
        (e.g. "species", "genus" or "family") and the values are the corresponding taxonomic names
    """
    lineage_data = gpf.ll(lineage_dump_path)
    lineage_dict = OrderedDict()
    for line in lineage_data:
        split_line = line.split("|")
        split_line = [n.strip() for n in split_line]
        entry_dict = {
            kraken_db_name + "_kraken_name": split_line[1],
            kraken_db_name + "_kraken_species": split_line[2],
            kraken_db_name + "_kraken_genus": split_line[3],
            kraken_db_name + "_kraken_family": split_line[4],
            kraken_db_name + "_kraken_order": split_line[5],
            kraken_db_name + "_kraken_class": split_line[6],
            kraken_db_name + "_kraken_phylum": split_line[7],
            kraken_db_name + "_kraken_kingdom": split_line[8],
            kraken_db_name + "_kraken_domain": split_line[9],
        }
        lineage_dict[split_line[0]] = entry_dict
    return lineage_dict


def get_kraken_and_lineage_dict(kraken_dict, lineage_dict, kraken_db_name):
    """
    Merges the kraken results with lineage information for the taxid numbers
    """
    kraken_and_lineage_dict = OrderedDict()
    for seq_name in kraken_dict:
        taxid = kraken_dict[seq_name]
        lineage_entry = {
            kraken_db_name + "_kraken_taxid": "0",
            kraken_db_name + "_kraken_name": None,
            kraken_db_name + "_kraken_species": None,
            kraken_db_name + "_kraken_genus": None,
            kraken_db_name + "_kraken_family": None,
            kraken_db_name + "_kraken_order": None,
            kraken_db_name + "_kraken_class": None,
            kraken_db_name + "_kraken_phylum": None,
            kraken_db_name + "_kraken_kingdom": None,
            kraken_db_name + "_kraken_domain": None,
        }
        if taxid in lineage_dict:
            lineage_entry = lineage_dict[taxid]
            lineage_entry[kraken_db_name + "_kraken_taxid"] = taxid
        else:
            if taxid != "0":
                sys.stderr.write("Taxid {} was not found in the lineages dump file\n".format(taxid))
        kraken_and_lineage_dict[seq_name] = lineage_entry

    return kraken_and_lineage_dict


def main(kraken_results_path, lineage_dump_path, kraken_db_name, out_path):
    kraken_dict = load_kraken_results(kraken_results_path)
    lineage_dict = load_lineage(lineage_dump_path, kraken_db_name)
    kraken_and_lineage_dict = get_kraken_and_lineage_dict(kraken_dict, lineage_dict, kraken_db_name)
    df = pd.DataFrame.from_dict(kraken_and_lineage_dict)
    df = df.transpose()
    df.index = df.index.rename("scaff")
    df.to_csv(out_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-v", "--version", action="version", version="1.0")
    parser.add_argument("kraken_results_path", help="Path to output file of a Kraken run", type=str)
    parser.add_argument(
        "lineage_dump_path",
        help="Path to an NCBI taxonomy rankedlineage.dmp file (downloaded from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz)",
        type=str,
    )
    parser.add_argument("kraken_db_name", help="Kraken database name", type=str, choices=["bacterial", "nt"])
    parser.add_argument("out_path", help="Path for output CSV file", type=str)
    args = parser.parse_args()
    main(args.kraken_results_path, args.lineage_dump_path, args.kraken_db_name, args.out_path)
