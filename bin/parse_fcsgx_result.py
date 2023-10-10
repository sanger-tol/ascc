#!/usr/bin/env python3
"""
Script for parsing the result files of FCS-GX, originally written by Eerik Aunin (ea10)
further minor modifications by Yumi Sims (yy5)
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
        sys.stderr.write("Could not find the directory with FCS-GX output files ({})\n".format(fcs_gx_reports_folder))
        sys.exit(1)

    fcs_gx_outfiles = os.listdir(fcs_gx_reports_folder)
    taxonomy_files = [n for n in fcs_gx_outfiles if n.endswith(".taxonomy.rpt")]
    if len(taxonomy_files) != 1:
        sys.stderr.write(
            "Error occurred when trying to find FCS-GX *.taxonomy.rpt file in {}\n".format(fcs_gx_reports_folder)
        )
        sys.exit(1)
    taxonomy_file = fcs_gx_reports_folder + "/" + taxonomy_files[0]

    report_files = [n for n in fcs_gx_outfiles if n.endswith(".fcs_gx_report.txt")]
    if len(report_files) != 1:
        sys.stderr.write(
            "Error occurred when trying to find FCS-GX *.fcs_gx_report.txt file in {}\n".format(fcs_gx_reports_folder)
        )
        sys.exit(1)
    report_file = fcs_gx_reports_folder + "/" + report_files[0]
    return taxonomy_file, report_file


def load_taxids_data(taxonomy_file):
    """
    Parses the *.taxonomy.rpt to find taxids that correspond to species names. Returns this as a dictionary
    """
    taxonomy_data = gpf.l(taxonomy_file)
    assert len(taxonomy_data) > 2
    taxonomy_data = taxonomy_data[2 : len(taxonomy_data)]
    collection_dict = OrderedDict()

    for line in taxonomy_data:
        add_new_entry = False
        multiple_divs_per_scaff = False
        split_line = line.split("\t")
        assert len(split_line) == 34
        scaff = split_line[0]
        if "~" in scaff:
            scaff = scaff.split("~")[0]
        tax_name_1 = split_line[5]
        tax_id_1 = split_line[6]
        div_1 = split_line[7]
        cvg_by_div_1 = split_line[8]
        if cvg_by_div_1 != "":
            cvg_by_div_1 = int(cvg_by_div_1)
        cvg_by_tax_1 = split_line[9]
        if cvg_by_tax_1 != "":
            cvg_by_tax_1 = int(cvg_by_tax_1)
        score_1 = split_line[10]
        if score_1 != "":
            score_1 = int(score_1)
        else:
            score_1 = 0
        if scaff in collection_dict:
            old_score = collection_dict[scaff]["fcs_gx_score"]
            if old_score < score_1:
                add_new_entry = True
                multiple_divs_per_scaff = True
        else:
            add_new_entry = True
        if add_new_entry is True:
            row_dict = dict()
            row_dict["fcs_gx_top_tax_name"] = tax_name_1
            row_dict["fcs_gx_top_taxid"] = tax_id_1
            row_dict["fcs_gx_div"] = div_1
            row_dict["fcs_gx_coverage_by_div"] = cvg_by_div_1
            row_dict["fcs_gx_coverage_by_tax"] = cvg_by_tax_1
            row_dict["fcs_gx_score"] = score_1
            row_dict["fcs_gx_multiple_divs_per_scaff"] = multiple_divs_per_scaff
            row_dict["fcs_gx_action"] = "NA"
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
    if os.path.isfile(rankedlineage_path) is False:
        sys.stderr.write(
            f"The NCBI rankedlineage.dmp file was not found at the expected location ({rankedlineage_path})\n"
        )
        sys.exit(1)
    rankedlineage_data = gpf.ll(rankedlineage_path)
    for line in rankedlineage_data:
        split_line = line.split("|")
        split_line = [n.strip() for n in split_line]
        assert len(split_line) == 11
        taxid = split_line[0]
        if taxid in taxids_list:
            current_lineage_dict = dict()
            for i in range(1, 10):
                current_lineage_dict[rankedlineage_col_names[i]] = split_line[i]
            lineages_dict[taxid] = current_lineage_dict
    return lineages_dict


def main(fcs_gx_reports_folder, ncbi_rankedlineage_path):
    taxonomy_file, report_file = get_fcs_gx_outfile_paths(fcs_gx_reports_folder)
    collection_dict = load_taxids_data(taxonomy_file)
    collection_dict = load_report_data(report_file, collection_dict)
    taxids_list = get_taxids_list(taxonomy_file)
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
        if row_top_taxid in lineages_dict:
            current_lineage_dict = lineages_dict[row_dict["fcs_gx_top_taxid"]]
            for i in range(1, 10):
                out_line += f",{current_lineage_dict[rankedlineage_col_names[i]]}"
                # current_lineage_dict[rankedlineage_col_names[i]] = split_line[i]
        else:
            for i in range(1, 10):
                out_line += ","
        print(out_line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "fcs_gx_reports_folder",
        type=str,
        help="Path to directory with FCS-GX output files (*.taxonomy.rpt and *.fcs_gx_report.txt)",
    )
    parser.add_argument("ncbi_rankedlineage_path", type=str, help="Path to the rankedlineage.dmp of NCBI taxonomy")
    parser.add_argument("-v", action="version", version="1.0")
    args = parser.parse_args()
    main(args.fcs_gx_reports_folder, args.ncbi_rankedlineage_path)
