#!/usr/bin/env python3

VERSION = "V1.0.0"

DESCRIPTION = """
-------------------------------------
            Autofilter
        Version = {VERSION}
-------------------------------------
Written by Eerik Aunin
Modified by Damon-Lee Pointon
-------------------------------------

Script for filtering the assembly to
remove putative contaminants based on
FCS-GX and Tiara results.
-------------------------------------

"""

from pathlib import Path
import general_purpose_functions as gpf
import os
import sys
import argparse
import textwrap


def parse_args():
    parser = argparse.ArgumentParser(
        prog="Autofilter",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(DESCRIPTION),
    )
    parser.add_argument("fasta", type=str, help="Path to the assembly FASTA file")
    parser.add_argument("-t", "--tiara", type=str, help="Path to the Tiara summary file")
    parser.add_argument("-s", "--fcsgx_summary", type=str, help="Path to the fcs-gx_summary.csv file")
    parser.add_argument(
        "-o",
        "--output_auto_filtered",
        type=str,
        help="Path to the assembly_autofiltered.fasta file",
        default="autofiltered.fasta",
    )
    #parser.add_argument(
    #    "-c", "--fcs_gx_and_tiara_summary", type=str, help="Path to the fcs-gx_and_tiara_combined_summary.csv file"
    #)
    parser.add_argument(
        "-r",
        "--rejected_seq",
        type=str,
        help="Path to the assembly_filtering_removed_sequences.txt file",
        default="assembly_filtering_removed_sequences.txt",
    )
    parser.add_argument("-i", "--taxid", type=int, help="NCBI taxonomy ID of the species")
    parser.add_argument(
        "-n", "--ncbi_rankedlineage_path", type=str, help="Path to the rankedlineage.dmp of NCBI taxonomy"
    )
    parser.add_argument(
        "--tiara_action_mode", type=str, choices=["warn", "remove"], default="warn", help="Action when Tiara detects a putative contaminant that is not reported as a contaminant by FCS-GX. The choices are 'warn' (print a warning) or 'remove' (remove this sequence from the assembly). Default: warn"
    )
    parser.add_argument("-v", "--version", action="version", version=VERSION)
    return parser.parse_args()


def get_domain_from_taxid(query_taxid, rankedlineage_path):
    """
    Input: 1) a taxID, 2) path to the NCBI rankedlineage.dmp file
    Output: domain classification corresponding to the taxID
    """
    domain = None
    query_taxid = str(query_taxid)
    rankedlineage_data = gpf.ll(rankedlineage_path)
    for line in rankedlineage_data:
        split_line = line.split("|")
        split_line = [n.strip() for n in split_line]
        assert len(split_line) == 11
        taxid = split_line[0]
        domain = split_line[9]
        if taxid == query_taxid:
            domain = split_line[9]
            if domain not in ("", "Archaea", "Bacteria", "Eukaryota", "Viruses"):
                sys.stderr.write(f"Unrecognised value for domain-level taxonomy: {domain}")
                sys.exit(1)
            break
    if domain is None:
        sys.stderr.write(
            "The domain for taxid ({}) was not found in the NCBI rankedlineage.dmp file ({})\n".format(
                query_taxid, rankedlineage_path
            )
        )
        sys.exit(1)
    return domain


def process_tiara_results(tiara_results_path, target_domain):
    """
    Input: 1) path to the main output file of Tiara, 2) the domain of the target species
    Output: dictionary where the keys are scaffold names and the values are the decontamination action based on Tiara results
        ('keep' or 'exclude')
    """
    tiara_action_dict = dict()

    allowed_classif_dict = dict()
    allowed_classif_dict[""] = ["archaea", "bacteria", "prokarya", "eukarya", "organelle", "unknown"]
    allowed_classif_dict["Archaea"] = ["archaea", "prokarya", "unknown"]
    allowed_classif_dict["Bacteria"] = ["bacteria", "prokarya", "unknown"]
    allowed_classif_dict["Eukaryota"] = ["eukarya", "organelle", "unknown"]
    allowed_classif_dict["Viruses"] = ["archaea", "bacteria", "prokarya", "eukarya", "organelle", "unknown"]
    allowed_classif_list = allowed_classif_dict[target_domain]

    tiara_output = gpf.ll(tiara_results_path)
    for counter, line in enumerate(tiara_output):
        if counter == 0:
            continue
        split_line = line.split()
        assert len(split_line) == 3
        tiara_class_fst_stage = split_line[1]
        assert tiara_class_fst_stage in ("archaea", "bacteria", "prokarya", "eukarya", "organelle", "unknown")
        tiara_action = "KEEP"
        if tiara_class_fst_stage not in allowed_classif_list:
            tiara_action = "EXCLUDE"
        scaff = split_line[0]
        tiara_action_dict[scaff] = tiara_action
    return tiara_action_dict


def get_fcs_gx_action_dict(fcs_gx_summary_path):
    """
    Input: path to FCS-GX summary CSV file (produced  by ascc_parse_fcsgx_results.py)
    Output: dictionary where the keys are scaffold names and the values are the FCS-GX action values
    """
    fcs_gx_action_dict = dict()
    fcs_gx_summary_data = gpf.ll(fcs_gx_summary_path)
    for counter, line in enumerate(fcs_gx_summary_data):
        if counter == 0:
            continue
        split_line = line.split(",")
        scaff = split_line[0]
        fcs_gx_action = split_line[8]
        fcs_gx_action_dict[scaff] = fcs_gx_action
    return fcs_gx_action_dict


def get_scaff_names(assembly_path):
    """
    Reads FASTA headers from a FASTA file and returns them as a list
    """
    scaffs = list()
    fasta_data = gpf.read_fasta_in_chunks(assembly_path)
    for fasta_tuple in fasta_data:
        scaffs.append(fasta_tuple[0])
    return scaffs


def filter_assembly(assembly_path, scaffs_to_exclude, filtered_assembly_path):
    """
    Filters a genome assembly FASTA file to remove sequences that are listed in the scaffs_to_exclude list
    """
    out_list = list()
    fasta_data = gpf.read_fasta_in_chunks(assembly_path)
    for header, seq in fasta_data:
        if header not in scaffs_to_exclude:
            out_list.append(">" + header)
            split_seq = gpf.split_with_fixed_row_length(seq, 80)
            out_list.extend(split_seq)
        else:
            sys.stderr.write(
                f"Excluding the sequence {header} from the filtered assembly ({filtered_assembly_path}), as it appears to be a contaminant based on FCS-GX and/or Tiara results\n"
            )
    gpf.export_list_as_line_break_separated_file(out_list, filtered_assembly_path)


def main():
    args = parse_args()
    if args.taxid == -1:
        sys.stderr.write(
            "The filtering of assembly based on FCS-GX and Tiara results requires a taxID but a valid taxID has not been provided (the provided taxID is -1, which is a placeholder value)\n"
        )

    assembly_path = args.fasta
    tiara_results_path = args.tiara
    fcs_gx_summary_path = args.fcsgx_summary
    filtered_assembly_path = args.output_auto_filtered
    #combined_summary = args.fcs_gx_and_tiara_summary
    excluded_seq_list_path = args.rejected_seq
    ncbi_rankedlist = args.ncbi_rankedlineage_path

    Path(f"./fasta/filtered").mkdir(parents=True, exist_ok=True)

    for i in [ncbi_rankedlist, tiara_results_path, fcs_gx_summary_path, assembly_path]:
        if not os.path.isfile(i):
            sys.stderr.write(f"{i} was not at the expected location\n")
            sys.exit(1)

    target_domain = get_domain_from_taxid(args.taxid, ncbi_rankedlist)
    tiara_action_dict = process_tiara_results(tiara_results_path, target_domain)
    fcs_gx_action_dict = get_fcs_gx_action_dict(fcs_gx_summary_path)

    combined_action_dict = dict()
    scaffs_to_exclude = list()
    scaffs = get_scaff_names(assembly_path)
    for scaff in scaffs:
        fcs_gx_action = "NA"
        tiara_action = "NA"
        if scaff in fcs_gx_action_dict:
            fcs_gx_action = fcs_gx_action_dict[scaff]
            combined_action_source = "FCS-GX"
        if scaff in tiara_action_dict:
            tiara_action = tiara_action_dict[scaff]
        combined_action = fcs_gx_action
        if fcs_gx_action == "NA" and tiara_action == "EXCLUDE":
            if args.tiara_action_mode == "remove":
                combined_action = "EXCLUDE"
                combined_action_source = "Tiara"
            elif args.tiara_action_mode == "warn":
                combined_action = "WARN"
                combined_action_source = "Tiara"
        if fcs_gx_action == "EXCLUDE" and tiara_action == "EXCLUDE":
            combined_action_source = "FCS-GX_and_Tiara"
        if combined_action == "EXCLUDE":
            scaffs_to_exclude.append(scaff)
        combined_action_dict[scaff] = {
            "fcs_gx_action": fcs_gx_action,
            "tiara_action": tiara_action,
            "combined_action": combined_action,
            "combined_action_source": combined_action_source,
        }
    filter_assembly(assembly_path, scaffs_to_exclude, filtered_assembly_path)
    gpf.export_list_as_line_break_separated_file(scaffs_to_exclude, excluded_seq_list_path)

    out_csv_list = list()
    out_csv_list.append("scaff,fcs_gx_action,tiara_action,combined_action")
    for scaff, scaff_properties in combined_action_dict.items():
        out_line = f"{scaff},{scaff_properties['fcs_gx_action']},{scaff_properties['tiara_action']},{scaff_properties['combined_action']},{scaff_properties['combined_action_source']}"
        out_csv_list.append(out_line)
    gpf.export_list_as_line_break_separated_file(out_csv_list, "ABNORMAL_CHECK.csv")


if __name__ == "__main__":
    main()