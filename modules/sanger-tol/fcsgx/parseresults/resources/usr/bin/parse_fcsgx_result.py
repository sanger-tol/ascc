#!/usr/bin/env python3
"""
Script for parsing the result files of FCS-GX, originally written by Eerik Aunin (ea10)
further refactoring/modifications by Yumi Sims (yy5)
Updated by Damon-Lee Pointon (dp24)
"""

import argparse
import os
import sys
from pathlib import Path


def file_to_generator(input_file: str):
    """
    Brought in from Eeriks GPF library

    Originally as gpf.l(input_file)
    """
    file_path = Path(input_file)

    if file_path.exists():
        with open(file_path) as f:
            for line in f:
                line = line.rstrip()
                yield line


def file_to_list(path: str) -> list:
    """
    Function for loading text file as a list and removing line breaks from line ends
    """
    try:
        return list(file_to_generator(path))
    except FileNotFoundError:
        sys.stderr.write("Error: file not found (" + path + ")\n")
        sys.exit(1)


def load_taxids_data(taxonomy_file: str) -> dict:
    """
    Parses the *.taxonomy.rpt to find taxids that correspond to species names. Returns this as a dictionary
    """
    taxonomy_data: list = file_to_list(taxonomy_file)[2:]
    assert taxonomy_data, "Taxonomy data is empty"

    collection_dict = dict()

    for line in taxonomy_data:
        split_line = line.split("\t")
        assert len(split_line) == 34

        scaff = split_line[0].split("~")[0]
        tax_name_1, tax_id_1, div_1, cvg_by_div_1, cvg_by_tax_1, score_1 = (
            split_line[5],
            split_line[6],
            split_line[7],
            int(split_line[8]) if split_line[8] else None,
            int(split_line[9]) if split_line[9] else None,
            int(split_line[10]) if split_line[10] else 0,
        )

        if scaff in collection_dict and int(collection_dict[scaff]["fcs_gx_score"] or 0) < score_1:
            row_dict = {
                "fcs_gx_top_tax_name": tax_name_1,
                "fcs_gx_top_taxid": tax_id_1,
                "fcs_gx_div": div_1,
                "fcs_gx_coverage_by_div": cvg_by_div_1,
                "fcs_gx_coverage_by_tax": cvg_by_tax_1,
                "fcs_gx_score": score_1,
                "fcs_gx_multiple_divs_per_scaff": True,
                "fcs_gx_action": "NA",
            }
            collection_dict[scaff] = row_dict
        elif scaff not in collection_dict:
            row_dict = {
                "fcs_gx_top_tax_name": tax_name_1,
                "fcs_gx_top_taxid": tax_id_1,
                "fcs_gx_div": div_1,
                "fcs_gx_coverage_by_div": cvg_by_div_1,
                "fcs_gx_coverage_by_tax": cvg_by_tax_1,
                "fcs_gx_score": score_1,
                "fcs_gx_multiple_divs_per_scaff": False,
                "fcs_gx_action": "NA",
            }
            collection_dict[scaff] = row_dict

    return collection_dict


def load_report_data(report_file: str, collection_dict: dict) -> dict:
    """
    Parses the *.fcs_gx_report.txt to add entries from the 'action' column to the collection of entries per scaffold that is stored in collection_dict
    """
    report_data = file_to_list(report_file)
    if len(report_data) > 2:
        report_data = report_data[2 : len(report_data)]
        for line in report_data:
            split_line = line.split("\t")
            assert len(split_line) == 8
            scaff: str = split_line[0]
            fcs_gx_action: str = split_line[4]
            collection_dict[scaff]["fcs_gx_action"] = fcs_gx_action
    return collection_dict


def get_taxids_list(fcs_gx_taxonomy_file_path: str) -> list:
    """
    Goes through FCS-GX taxonomy output file and returns a list of unique taxIDs found in the file
    """
    if not os.path.isfile(fcs_gx_taxonomy_file_path):
        sys.stderr.write(
            f"The FCS-GX taxonomy file was not found at the expected location ({fcs_gx_taxonomy_file_path})\n"
        )
        sys.exit(1)
    taxids_list = list()
    for line in file_to_list(fcs_gx_taxonomy_file_path):
        if not line.startswith("#"):
            split_line = line.split("\t")
            assert len(split_line) == 34
            taxid = split_line[6]
            if taxid not in taxids_list:
                taxids_list.append(taxid)
    return taxids_list


def get_lineages_by_taxid(taxids_list: list, rankedlineage_path: str) -> dict:
    """
    Takes a list of taxIDs and the path to the NCBI rankedlineage.dmp file as the input. Returns the lineage corresponding to each taxID
    """
    lineages_dict = dict()
    rankedlineage_col_names = (
        "taxid",
        "fcs_gx_name",
        "fcs_gx_species",
        "fcs_gx_genus",
        "fcs_gx_family",
        "fcs_gx_order",
        "fcs_gx_class",
        "fcs_gx_phylum",
        "fcs_gx_kingdom",
        "fcs_gx_domain",
    )

    for line in file_to_list(rankedlineage_path):
        split_line: list = line.split("|")
        split_line: list = [n.strip() for n in split_line]
        assert len(split_line) >= 11, (
            f"Expected at least 11 columns in rankedlineage.dmp, got {len(split_line)}"
        )  # This should now handle both new and old formats (as of April 2025)
        taxid = split_line[0]
        if taxid in taxids_list:
            current_lineage_dict = dict()
            for i in range(1, 10):
                current_lineage_dict[rankedlineage_col_names[i]] = split_line[i]
            lineages_dict[taxid] = current_lineage_dict
    return lineages_dict


def main(taxonomy_report: str, fcs_report: str, ncbi_rankedlineage_path: str) -> None:
    collection_dict = load_taxids_data(taxonomy_report)
    collection_dict = load_report_data(fcs_report, collection_dict)
    taxids_list = get_taxids_list(taxonomy_report)
    lineages_dict = get_lineages_by_taxid(taxids_list, ncbi_rankedlineage_path)
    rankedlineage_col_names = (
        "taxid",
        "fcs_gx_name",
        "fcs_gx_species",
        "fcs_gx_genus",
        "fcs_gx_family",
        "fcs_gx_order",
        "fcs_gx_class",
        "fcs_gx_phylum",
        "fcs_gx_kingdom",
        "fcs_gx_domain",
    )

    out_header = "scaff,fcs_gx_top_tax_name,fcs_gx_top_taxid,fcs_gx_div,fcs_gx_coverage_by_div,fcs_gx_coverage_by_tax,fcs_gx_score,fcs_gx_multiple_divs_per_scaff,fcs_gx_action"
    out_header += "," + ",".join(rankedlineage_col_names)
    print(out_header)
    for scaff, row_dict in collection_dict.items():
        out_line = [
            scaff,
            row_dict["fcs_gx_top_tax_name"],
            row_dict["fcs_gx_top_taxid"],
            row_dict["fcs_gx_div"],
            row_dict["fcs_gx_coverage_by_div"],
            row_dict["fcs_gx_coverage_by_tax"],
            row_dict["fcs_gx_score"],
            row_dict["fcs_gx_multiple_divs_per_scaff"],
            row_dict["fcs_gx_action"],
        ]
        row_top_taxid: str = row_dict["fcs_gx_top_taxid"]
        if row_top_taxid in lineages_dict:
            current_lineage_dict: dict = lineages_dict[row_dict["fcs_gx_top_taxid"]]
            for i in range(1, 10):
                out_line.append(current_lineage_dict[rankedlineage_col_names[i]])
        else:
            for _ in range(1, 10):
                out_line.append("")
        print(*out_line, sep=",")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "taxonomy_report",
        type=str,
        help="Path to directory with FCS-GX output files *.taxonomy.rpt",
    )
    parser.add_argument(
        "fcs_report",
        type=str,
        help="Path to directory with FCS-GX output files *.fcs_gx_report.txt",
    )
    parser.add_argument(
        "ncbi_rankedlineage_path",
        type=str,
        help="Path to the rankedlineage.dmp of NCBI taxonomy",
    )
    parser.add_argument("--version", action="version", version="1.1.0")
    args = parser.parse_args()

    main(args.taxonomy_report, args.fcs_report, args.ncbi_rankedlineage_path)
