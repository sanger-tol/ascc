#!/usr/bin/env python3
"""
Script for parsing the result files of FCS-GX, originally written by Eerik Aunin (ea10)
further refactoring/modifications by Yumi Sims (yy5)
Updated by Damon-Lee Pointon (dp24)
"""

import general_purpose_functions as gpf
from collections import OrderedDict
import os
import sys
import argparse


def get_fcs_gx_outfile_paths(fcs_gx_reports_folder):
    """
    Locates the FCS-GX output files in the FCS-GX output directory
    """
    if os.path.isdir(fcs_gx_reports_folder) is False:
        sys.stderr.write(
            "Could not find the directory with FCS-GX output files ({})\n".format(
                fcs_gx_reports_folder
            )
        )
        sys.exit(1)

    fcs_gx_outfiles = os.listdir(fcs_gx_reports_folder)
    taxonomy_files = [n for n in fcs_gx_outfiles if n.endswith(".taxonomy.rpt")]
    if len(taxonomy_files) != 1:
        sys.stderr.write(
            "Error occurred when trying to find FCS-GX *.taxonomy.rpt file in {}\n".format(
                fcs_gx_reports_folder
            )
        )
        sys.exit(1)
    taxonomy_file = fcs_gx_reports_folder + "/" + taxonomy_files[0]

    report_files = [n for n in fcs_gx_outfiles if n.endswith(".fcs_gx_report.txt")]
    if len(report_files) != 1:
        sys.stderr.write(
            "Error occurred when trying to find FCS-GX *.fcs_gx_report.txt file in {}\n".format(
                fcs_gx_reports_folder
            )
        )
        sys.exit(1)
    report_file = fcs_gx_reports_folder + "/" + report_files[0]
    return taxonomy_file, report_file


def load_taxids_data(taxonomy_file):
    """
    Parses the *.taxonomy.rpt to find taxids that correspond to species names. Returns this as a dictionary
    """
    taxonomy_data = gpf.l(taxonomy_file)[2:]
    assert taxonomy_data, "Taxonomy data is empty"

    collection_dict = OrderedDict()

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

        if (
            scaff in collection_dict
            and collection_dict[scaff]["fcs_gx_score"] < score_1
        ):
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


def load_report_data(report_file, collection_dict):
    """
    Parses the *.fcs_gx_report.txt to add entries from the 'action' column to the collection of entries per scaffold that is stored in collection_dict
    """
    report_data = gpf.l(report_file)
    if len(report_data) > 2:
        report_data = report_data[2 : len(report_data)]
        for line in report_data:
            split_line = line.split("\t")
            assert len(split_line) == 8
            scaff = split_line[0]
            fcs_gx_action = split_line[4]
            collection_dict[scaff]["fcs_gx_action"] = fcs_gx_action
    return collection_dict


def get_taxids_list(fcs_gx_taxonomy_file_path):
    """
    Goes through FCS-GX taxonomy output file and returns a list of unique taxIDs found in the file
    """
    if os.path.isfile(fcs_gx_taxonomy_file_path) is False:
        sys.stderr.write(
            f"The FCS-GX taxonomy file was not found at the expected location ({fcs_gx_taxonomy_file_path})\n"
        )
        sys.exit(1)
    taxids_list = list()
    fcs_gx_taxonomy_data = gpf.ll(fcs_gx_taxonomy_file_path)
    for line in fcs_gx_taxonomy_data:
        if line.startswith("#") is False:
            split_line = line.split("\t")
            assert len(split_line) == 34
            taxid = split_line[6]
            if taxid not in taxids_list:
                taxids_list.append(taxid)
    return taxids_list


def get_lineages_by_taxid(taxids_list, rankedlineage_path):
    """
    Takes a list of taxIDs and the path to the NCBI rankedlineage.dmp file as the input. Returns the lineage corresponding to each taxID
    """
    lineages_dict = dict()
    # The columns in rankedlineage.dmp are: tax_id, tax_name, species, genus, family, order, class, phylum, kingdom, realm, domain
    if os.path.isfile(rankedlineage_path) is False:
        sys.stderr.write(
            f"The NCBI rankedlineage.dmp file was not found at the expected location ({rankedlineage_path})\n"
        )
        sys.exit(1)

    # Log the number of taxids we're looking for
    sys.stderr.write(f"Looking for {len(taxids_list)} taxids in rankedlineage.dmp\n")
    if len(taxids_list) > 0:
        sys.stderr.write(f"First few taxids: {', '.join(taxids_list[:5])}\n")

    line_count = 0
    found_taxids = 0

    for line in gpf.ll(rankedlineage_path):
        line_count += 1
        split_line = line.split("|")
        split_line = [n.strip() for n in split_line]
        if len(split_line) < 11:
            sys.stderr.write(
                f"Warning: Not enough columns in rankedlineage.dmp line, got {len(split_line)}, expected at least 11. Skipping line.\n"
            )
            continue

        taxid = split_line[0]
        if taxid in taxids_list:
            found_taxids += 1
            current_lineage_dict = dict()
            current_lineage_dict["taxid"] = split_line[0]
            current_lineage_dict["fcs_gx_name"] = (
                split_line[1] if len(split_line) > 1 else ""
            )
            current_lineage_dict["fcs_gx_species"] = (
                split_line[2] if len(split_line) > 2 else ""
            )
            current_lineage_dict["fcs_gx_genus"] = (
                split_line[3] if len(split_line) > 3 else ""
            )
            current_lineage_dict["fcs_gx_family"] = (
                split_line[4] if len(split_line) > 4 else ""
            )
            current_lineage_dict["fcs_gx_order"] = (
                split_line[5] if len(split_line) > 5 else ""
            )
            current_lineage_dict["fcs_gx_class"] = (
                split_line[6] if len(split_line) > 6 else ""
            )
            current_lineage_dict["fcs_gx_phylum"] = (
                split_line[7] if len(split_line) > 7 else ""
            )
            current_lineage_dict["fcs_gx_kingdom"] = (
                split_line[8] if len(split_line) > 8 else ""
            )
            current_lineage_dict["fcs_gx_realm"] = (
                split_line[9] if len(split_line) > 9 else ""
            )
            current_lineage_dict["fcs_gx_domain"] = (
                split_line[10] if len(split_line) > 10 else ""
            )

            lineages_dict[taxid] = current_lineage_dict

    # Log how many taxids we found and total lines processed
    sys.stderr.write(f"Processed {line_count} lines from rankedlineage.dmp\n")
    sys.stderr.write(
        f"Found {found_taxids} out of {len(taxids_list)} taxids in rankedlineage.dmp\n"
    )

    return lineages_dict


def main(fcs_gx_reports_folder, ncbi_rankedlineage_path):
    taxonomy_file, report_file = get_fcs_gx_outfile_paths(fcs_gx_reports_folder)
    collection_dict = load_taxids_data(taxonomy_file)
    collection_dict = load_report_data(report_file, collection_dict)
    taxids_list = get_taxids_list(taxonomy_file)
    lineages_dict = get_lineages_by_taxid(taxids_list, ncbi_rankedlineage_path)
    # The columns in rankedlineage.dmp are: tax_id, tax_name, species, genus, family, order, class, phylum, kingdom, realm, domain
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
        "fcs_gx_realm",
        "fcs_gx_domain",
    )

    # Log the number of scaffolds and taxids
    sys.stderr.write(f"Processing {len(collection_dict)} scaffolds\n")
    sys.stderr.write(f"Found {len(lineages_dict)} taxids with lineage information\n")

    # Check if any taxids have empty strings
    empty_taxids = [taxid for taxid in taxids_list if not taxid]
    if empty_taxids:
        sys.stderr.write(f"Warning: Found {len(empty_taxids)} empty taxids\n")

    out_header = "scaff,fcs_gx_top_tax_name,fcs_gx_top_taxid,fcs_gx_div,fcs_gx_coverage_by_div,fcs_gx_coverage_by_tax,fcs_gx_score,fcs_gx_multiple_divs_per_scaff,fcs_gx_action"
    out_header += "," + ",".join(rankedlineage_col_names)
    print(out_header)

    # Count how many scaffolds have lineage information
    scaffolds_with_lineage = 0

    for scaff, row_dict in collection_dict.items():
        out_line = "{},{},{},{},{},{},{},{},{}".format(
            scaff,
            row_dict["fcs_gx_top_tax_name"],
            row_dict["fcs_gx_top_taxid"],
            row_dict["fcs_gx_div"],
            row_dict["fcs_gx_coverage_by_div"],
            row_dict["fcs_gx_coverage_by_tax"],
            row_dict["fcs_gx_score"],
            row_dict["fcs_gx_multiple_divs_per_scaff"],
            row_dict["fcs_gx_action"],
        )
        row_top_taxid = row_dict["fcs_gx_top_taxid"]

        # Log if the taxid is empty
        if not row_top_taxid:
            sys.stderr.write(f"Warning: Empty taxid for scaffold {scaff}\n")

        if row_top_taxid and row_top_taxid in lineages_dict:
            scaffolds_with_lineage += 1
            current_lineage_dict = lineages_dict[row_dict["fcs_gx_top_taxid"]]
            # Iterate through the rankedlineage_col_names to get values from the dictionary
            for col_name in rankedlineage_col_names:
                out_line += f",{current_lineage_dict.get(col_name, '')}"  # Use .get to handle missing realms/domains gracefully
        else:
            # Add empty columns for all rankedlineage_col_names if taxid not found
            for _ in rankedlineage_col_names:
                out_line += ","
        print(out_line)

    # Log how many scaffolds have lineage information
    sys.stderr.write(
        f"Found lineage information for {scaffolds_with_lineage} out of {len(collection_dict)} scaffolds\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "fcs_gx_reports_folder",
        type=str,
        help="Path to directory with FCS-GX output files (*.taxonomy.rpt and *.fcs_gx_report.txt)",
    )
    parser.add_argument(
        "ncbi_rankedlineage_path",
        type=str,
        help="Path to the rankedlineage.dmp of NCBI taxonomy",
    )
    parser.add_argument("-v", action="version", version="1.0.1")
    args = parser.parse_args()
    main(args.fcs_gx_reports_folder, args.ncbi_rankedlineage_path)
