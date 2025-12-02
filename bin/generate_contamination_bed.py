#!/usr/bin/env python

# Adapted from ascc_nf_core_to_bed.py by James Torrance (jt8@sanger.ac.uk)

from typing import Dict, List, Tuple
import os
import sys
import re
import collections
import csv
import gzip
from Bio import SeqIO
import argparse

"""Generate BED file from decon results"""

def parse_args(args=None):
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--decon_tolid_type_dir",
        type=str, 
        help=".",
    )
    parser.add_argument(
        "--assembly_type",
        type=str, 
        help="."
    )
    parser.add_argument(
		"--assembly_path", 
		type=str, 
		help=""
	)
    parser.add_argument(
        "--is_organelle", 
        type=bool, 
        help="Is an organelle", 
    )
    parser.add_argument(
        "--no_fcs_gx", 
        type=bool, 
        help="Do not use FCS-GX"
    )
    parser.add_argument(
        "--no_barcodes",
        help="No multiplexing barcodes expected",
    )
    parser.add_argument(
        "--longread_paths", 
        type=str, 
        nargs="?",
        help="", 
    )
    parser.add_argument("--version", action="version", version="%(prog)s 1.0")
    
    return parser.parse_args(args)

def main():
    args = parse_args()
    
    merged_filter_file, fcs_gx_file, fcs_adaptor_file, trim_ns_file, barcode_files = get_decon_bed_files(decon_tolid_type_dir, is_organelle)

    # BASE BED FILE NAMES OFF ASSEMBLY
    main_bed_output_file = f"{re.sub(".fa(sta)?.gz$", "", args.assembly_path)}.contamination.bed"
    tiara_bed_output_file = f"{re.sub(".fa(sta)?.gz$", "", args.assembly_path)}.tiara.bed"

    main_bed_output_handle = open(main_bed_output_file, "w")
    tiara_bed_output_handle = open(tiara_bed_output_file, "w")

    # Parse FASTA assembly
    length_for_scaffold, scaffold_count = get_scaffold_for_fasta(args.assembly_path)
    lengths_removed = []

    # Parse FCS-GX file
    fcs_gx_taxonomy_for_scaffold = {}
    if not args.no_fcs_gx and not args.is_organelle:
        fcs_gx_taxonomy_for_scaffold = parse_fcs_gx_file(fcs_gx_file)

    sequences_removed_by_reason = collections.defaultdict(list)
    scaffolds_removed = 0

    # Parse FCS-Adaptor file
    main_bed_output_handle.write("# FCS-Adaptor\n")
    fcs_sequences_removed = parse_fcs_adapator_file(fcs_adaptor_file, length_for_scaffold)
    for key, val in fcs_sequences_removed.items():
        main_bed_output_handle.write(key)
        sequences_removed_by_reason["FCS-Adaptor"].append(val)

    long_organelle_name = {
        "mito": "Mitochondrion",
        "plastid": "Plastid",
    }

    # Parse organelle csv files
    if not args.is_organelle:
        for organelle in ["mito", "plastid"]:
            organelle_recommendation_dir = search_directory_for_file_pattern(
                args.decon_tolid_type_dir,
                f"organelle_contamination_recommendations.*_{organelle.upper()}$",
                no_files_okay=True,
            )
            if organelle_recommendation_dir is not None and os.path.isdir(
                organelle_recommendation_dir
            ):
                organelle_file = search_directory_for_file_pattern(
                    organelle_recommendation_dir,
                    f"{organelle.upper()}.*.bed$",
                    no_files_okay=True,
                )
                if organelle_file is not None:
                    main_bed_output_handle.write(f"# {long_organelle_name[organelle]}\n")
                    organelle_sequences_removed, scaffolds_removed_for_organelle = parse_organelle_csv(organelle_file)
                    for key, val in organelle_sequences_removed.items():
                        main_bed_output_handle.write(key)
                        sequences_removed_by_reason[long_organelle_name[organelle]].append(val)
                        lengths_removed.append(val)
                    scaffolds_removed += scaffolds_removed_for_organelle

    # Parse Barcode files
    if not args.no_barcodes:
        main_bed_output_handle.write("# Barcodes" + "\n")
        barcode_sequences_removed = parse_barcode_files(barcode_files)
        for key, val in barcode_sequences_removed.items():
            main_bed_output_handle.write(key)
            sequences_removed_by_reason["Barcodes"].append(val)
            lengths_removed.append(val)

    # Parse Trim Ns file
    main_bed_output_handle.write("# Trim Ns\n")
    trim_ns_sequences_removed = parse_trim_ns_file(trim_ns_file)
    for key, val in trim_ns_sequences_removed.items():
        main_bed_output_handle.write(key)
        sequences_removed_by_reason["Trim Ns"].append(val)
        lengths_removed.append(val)

    fcs_gx_taxonomy_removed = collections.defaultdict(list)
    fcs_ambiguity_text = ""

    if not args.no_fcs_gx and not args.is_organelle:
        main_bed_output_handle.write("# Merged ASCC call\n")
        with open(merged_filter_file) as merged_filter_handle:
            for line in merged_filter_handle:
                fields = re.split(",", line.rstrip())
                if fields[3] == "EXCLUDE":
                    scaffold = fields[0]
                    start = 0
                    if fields[0] in length_for_scaffold:
                        end = length_for_scaffold[fields[0]]
                    else:
                        exit(f"Cannot find scaffold {scaffold} in FASTA assembly")
                    treatment = "REMOVE"
                    if fields[4] in ("FCS-GX_and_Tiara", "FCS-GX"):
                        main_bed_output_handle.write(
                            "\t".join([scaffold, str(start), str(end), treatment])
                            + "\n"
                        )
                        lengths_removed.append(end)

                        scaffolds_removed += 1

                        if scaffold in fcs_gx_taxonomy_for_scaffold:
                            fcs_gx_taxonomy_removed[
                                fcs_gx_taxonomy_for_scaffold[scaffold]
                            ].append(length_for_scaffold[fields[0]])
                            # TODO Put in an exception for scaffolds where FCS-GX cannot make a determination
                        else:
                            exit(f"Cannot find {scaffold} in FCS-GX file")
                    elif fields[4] == "Tiara":
                        tiara_bed_output_handle.write(
                            "\t".join([scaffold, str(start), str(end), treatment])
                            + "\n"
                        )
                # Cases called as REVIEW or INFO by FCS-GX should be manually reviewed
                if fields[1] in ("REVIEW", "INFO"):
                    fcs_ambiguity_text += f"FCS-GX assigns scaffold {scaffold} ambiguous verdict {fields[1]} which requires manual review\n"

    main_bed_output_handle.close()
    tiara_bed_output_handle.close()

    contamination_report, is_abnormal, abnormal_details = build_contamination_report(args.assembly_type, length_for_scaffold, lengths_removed, scaffolds_removed, scaffold_count, fcs_gx_taxonomy_removed, sequences_removed_by_reason, fcs_ambiguity_text)

    return contamination_report, is_abnormal, abnormal_details


def build_contamination_report(assembly_type: str, length_for_scaffold, lengths_removed, scaffolds_removed, scaffold_count, fcs_gx_taxonomy_removed, sequences_removed_by_reason, fcs_ambiguity_text: str   ):
    """Build contamination report for JIRA ticket, and determine if abnormal_contamination_report label should be added"""

    report_is_abnormal = False
    abnormal_details = ""

    total_assembly_length = sum(length_for_scaffold.values())
    value_for_parameter = {
        "TOTAL_LENGTH_REMOVED": sum(lengths_removed),
        "PERCENTAGE_LENGTH_REMOVED": 100 * sum(lengths_removed) / total_assembly_length,
        "LARGEST_SCAFFOLD_REMOVED": max(lengths_removed, default=0),
        "SCAFFOLDS_REMOVED": scaffolds_removed,
        "PERCENTAGE_SCAFFOLDS_REMOVED": 100 * scaffolds_removed / scaffold_count,
    }

    jira_report = f"Contamination report for assembly labelled {assembly_type}\n"
    jira_report += f"Total length of scaffolds removed: {value_for_parameter['TOTAL_LENGTH_REMOVED']:,} ({value_for_parameter['PERCENTAGE_LENGTH_REMOVED']:.1f} %)\n"
    jira_report += f"Scaffolds removed: {value_for_parameter['SCAFFOLDS_REMOVED']} ({value_for_parameter['PERCENTAGE_SCAFFOLDS_REMOVED']:.1f} %)\n"
    jira_report += f"Largest scaffold removed: ({value_for_parameter['LARGEST_SCAFFOLD_REMOVED']:,})\n"
    jira_report += "\nFCS-GX contaminant species (number of scaffolds; total length of scaffolds): \n"
    fcs_gx_taxonomy_removed = dict(
        sorted(
            fcs_gx_taxonomy_removed.items(), key=lambda item: len(item[1]), reverse=True
        )
    )
    for taxonomy in fcs_gx_taxonomy_removed:
        jira_report += f"\t{taxonomy} ({len(fcs_gx_taxonomy_removed[taxonomy])}; {sum(fcs_gx_taxonomy_removed[taxonomy]):,})\n"

    alarm_threshold_for_parameter = {
        "TOTAL_LENGTH_REMOVED": 1e7,
        "PERCENTAGE_LENGTH_REMOVED": 3,
        "LARGEST_SCAFFOLD_REMOVED": 2e6,
        "PERCENTAGE_SCAFFOLDS_REMOVED": 10,  # ???
    }

    alarm_threshold_for_parameter_by_reason = {
        #'Mitochondrion': {'SCAFFOLDS_REMOVED': 10},
        #'Plastid': {'SCAFFOLDS_REMOVED': 10},
    }

    value_for_parameter_by_reason = {}

    for reason in ("Mitochondrion", "Plastid", "FCS-Adaptor", "Barcodes", "Trim Ns"):
        alarm_by_reason_description = ""
        if reason in sequences_removed_by_reason:
            value_for_parameter_by_reason[reason] = {}
            jira_report += f"\n{reason} ({len(sequences_removed_by_reason[reason])}; {sum(sequences_removed_by_reason[reason]):,})\n"
            if reason in ("Mitochondrion", "Plastid"):
                alarm_threshold_for_parameter["TOTAL_LENGTH_REMOVED"] += sum(
                    sequences_removed_by_reason[reason]
                )
                alarm_threshold_for_parameter["PERCENTAGE_LENGTH_REMOVED"] += (
                    100
                    * sum(sequences_removed_by_reason[reason])
                    / total_assembly_length
                )
                alarm_threshold_for_parameter["PERCENTAGE_SCAFFOLDS_REMOVED"] += (
                    100 * len(sequences_removed_by_reason[reason]) / scaffold_count
                )
            value_for_parameter_by_reason[reason]["SCAFFOLDS_REMOVED"] = len(
                sequences_removed_by_reason[reason]
            )
            for parameter in value_for_parameter_by_reason[reason]:
                if (
                    reason in alarm_threshold_for_parameter_by_reason
                    and parameter in alarm_threshold_for_parameter_by_reason[reason]
                    and value_for_parameter_by_reason[reason][parameter]
                    >= alarm_threshold_for_parameter_by_reason[reason][parameter]
                ):
                    alarm_by_reason_description += f" in assembly {assembly_type} {parameter} is {value_for_parameter_by_reason[reason][parameter]} which is above the alarm threshold {alarm_threshold_for_parameter_by_reason[reason][parameter]}."
                    print("abnormal_contamination_report triggered")
                    report_is_abnormal = True
                    print(alarm_by_reason_description)

    for parameter in value_for_parameter:
        alarm_description = ""
        print(f"{parameter}: {value_for_parameter[parameter]}")
        if (
            parameter in alarm_threshold_for_parameter
            and value_for_parameter[parameter]
            >= alarm_threshold_for_parameter[parameter]
        ):
            alarm_description += f" in assembly {assembly_type} {parameter} is {value_for_parameter[parameter]} which is above the alarm threshold {alarm_threshold_for_parameter[parameter]}."
            print("abnormal_contamination_report triggered")
            report_is_abnormal = True
            print(alarm_description)
        if alarm_description != "":
            alarm_description = "Abnormal contamination report:" + alarm_description
            abnormal_details += alarm_description
            jira_report += f"\n{alarm_description}\n"

    if fcs_ambiguity_text != "":
        jira_report += f"\n{fcs_ambiguity_text}\n"
        report_is_abnormal = True
        abnormal_details += fcs_ambiguity_text

    return jira_report, report_is_abnormal, abnormal_details

def parse_fcs_gx_file(fcs_gx_file: str) -> Dict[str, str]:
    fcs_gx_taxonomy_for_scaffold = {}
    with open(fcs_gx_file) as fcs_gx_handle:
        next(fcs_gx_handle)  # Skip header
        for line in fcs_gx_handle:
            fields = re.split(",", line.rstrip())
            if len(fields) < 4:
                exit(f"FCS-GX - Not enough fields in line {line}")
                sys
            scaffold = fields[0]
            taxonomy = fields[1]
            fcs_gx_div = fields[3]
            if re.search(":", fcs_gx_div):
                (fcs_gx_div_short, fcs_gx_div_long) = re.split(":", fcs_gx_div)
                taxonomy = f"{taxonomy}, {fcs_gx_div_long}"
            fcs_gx_taxonomy_for_scaffold[scaffold] = taxonomy

    return fcs_gx_taxonomy_for_scaffold

def parse_fcs_adapator_file(fcs_adaptor_file: str, length_for_scaffold) -> Dict[str, int]:
    sequences_removed = {}
    with open(fcs_adaptor_file) as fcs_adaptor_handle:
        for line in fcs_adaptor_handle:
            fields = re.split(r"\s+", line.rstrip())
            if not re.match("^#", line):
                accession = fields[0]
                length = fields[1]
                action = fields[2]
                ranges = fields[3]
                start = None
                end = None
                treatment = None
                if action == "ACTION_TRIM":
                    range_list = re.split(",", ranges)
                    for range in range_list:
                        (start, end) = re.split(r"\.\.", range)
                        start = int(start)
                        end = int(end)
                        start -= 1
                        treatment = "CONTAMINANT"
                # If only a short (< ~200 bp) sequence is left after trimming, FCS-Adaptor excludes it
                elif action == "ACTION_EXCLUDE":
                    start = 0
                    if accession in length_for_scaffold:
                        end = length_for_scaffold[accession]
                    else:
                        exit(f"Cannot find scaffold {accession} in FASTA assembly")
                    treatment = "REMOVE"
                else:
                    exit(f"Action not recognised: {action}")

                main_output = "\t".join([accession, str(start), str(end), treatment]) + "\n"
                sequences_removed[main_output] = end - start
        return sequences_removed

def parse_barcode_files(barcode_files: List[str]) -> Dict[str, int]:
    sequences_removed = {}        
    for barcode_file in barcode_files:
        with open(barcode_file) as barcode_handle:
            for line in barcode_handle:
                fields = re.split(r"\s+", line.rstrip())
                if fields[0] == "CONTAMINANT":
                    scaffold = fields[1]
                    start = int(fields[2]) - 1
                    end = int(fields[3])
                    treatment = "CONTAMINANT"
                    main_output = "\t".join([scaffold, str(start), str(end), treatment])+ "\n"
                    sequences_removed[main_output] = end - start
    return sequences_removed

def parse_trim_ns_file(trim_ns_file: str) -> Dict[str, int]:
    sequences_removed = {}
    with open(trim_ns_file) as trim_ns_handle:
        for line in trim_ns_handle:
            fields = re.split(r"\s+", line.rstrip())
            if fields[0] in (
                "TRIM",
                "REVCLIP",
                "FWDCLIP",
                "CLIP",
            ):  # You may just want one
                scaffold = fields[1]
                start = int(fields[2]) - 1
                end = fields[3]
                treatment = "TRIM"
                main_output ="\t".join([scaffold, str(start), str(end), treatment]) + "\n"
                sequences_removed[main_output] = end - start
    return sequences_removed

def parse_organelle_csv(organelle_file: str):
    sequences_removed = {}
    scaffolds_removed = 0
    with open(organelle_file) as organelle_handle:
        bed_csv_reader = csv.reader(
            organelle_handle, delimiter="\t"
        )
        for field_set in bed_csv_reader:
            scaffold = field_set[0]
            start = int(field_set[1])
            end = int(field_set[2])
            (length, percentage) = re.split(",", field_set[3])
            percentage = float(percentage)
            # Discard anything over 95% coverage, or over 90% coverage if it's below 50 kb
            if (percentage > 90 and end < 50000) or (
                percentage > 95
            ):
                treatment = "REMOVE"
                main_output = "\t".join([scaffold, str(start), str(end), treatment])+ "\n"
                sequences_removed[main_output] = end
    return sequences_removed, scaffolds_removed

def get_decon_bed_files(decon_tolid_type_dir: str, is_organelle: bool, no_fcs_gx: bool, no_barcodes: bool, longread_paths: List[str]) -> Dict[str, str]: 
    if not no_fcs_gx and not is_organelle:
        merged_filter_file = search_directory_for_file_pattern(
            decon_tolid_type_dir + "/autofilter/", "_ABNORMAL_CHECK.csv$"
        )
        fcs_gx_file = search_directory_for_file_pattern(
            decon_tolid_type_dir + "/fcsgx_data/", "_parsed_fcsgx.csv$"
        )
    else:
        merged_filter_file = None
        fcs_gx_file = None
        
    fcs_adaptor_file = search_directory_for_file_pattern(
        decon_tolid_type_dir + "/fcs_adaptor/", "_euk.fcs_adaptor_report.txt$"
    )
    trim_ns_file = search_directory_for_file_pattern(
        decon_tolid_type_dir + "/trailingns/", "_trim_Ns$"
    )
    if trim_ns_file is None:
        print("Could not find trim Ns file")
        sys.exit(1)

    barcode_files = search_directory_for_file_pattern(
        decon_tolid_type_dir + "/filter_barcode/", "_filtered.txt$", multiple_results=True
    )

    if len(barcode_files) == 0 and not no_barcodes:
        # Were there barcodes in the input file?
        barcodes_found = True
        for fasta_pacbio_read_dir in longread_paths:
            pacbio_reads_file = None
            if os.path.isdir(fasta_pacbio_read_dir):
                pacbio_reads_files = search_directory_for_file_pattern(
                    fasta_pacbio_read_dir, ".fa(sta)?.gz$", multiple_results=True
                )
                # Quit if there was a barcode in the input filename, but otherwise continue
                for pacbio_reads_file in pacbio_reads_files:
                    if re.search(r"bc\d{4}", pacbio_reads_file):
                        barcodes_found = False
            else:
                print(
                    f"Cannot find FASTA subdir of pacbio read dir: nothing at {fasta_pacbio_read_dir}"
                )
                sys.exit(1)

        if not barcodes_found:
            print(f"Could not find barcode file in {decon_tolid_type_dir}")
            sys.exit(1)

    return merged_filter_file, fcs_gx_file, fcs_adaptor_file, trim_ns_file, barcode_files

def search_directory_for_file_pattern(directory, pattern, multiple_results=False, no_files_okay = False):
    if os.path.isdir(directory):
        match_files = []
        for candidate_file in os.scandir(directory):
            if re.search(pattern, candidate_file.name):
                match_files.append(directory + '/' + candidate_file.name)
        if multiple_results:
            return match_files
        elif len(match_files) == 1:
            return(match_files[0])
        else:
            if len(match_files) == 0 and no_files_okay:
                return None
            else:
                print(f'Could not find unique file for pattern {pattern}: candidates include {match_files}')
                sys.exit(1)
    else:
        print(f'Not a directory: {directory}')
        sys.exit(1)


def get_scaffold_for_fasta(assembly_file: str) -> Tuple[Dict[str, int], int]:
    """Get list of scaffolds from fasta file"""

    fasta_input_handle = gzip.open(assembly_file, "rt")
    length_for_scaffold = {}
    scaffold_count = 0
    for record in SeqIO.parse(fasta_input_handle, "fasta"):
        length_for_scaffold[record.id] = len(record.seq)
        scaffold_count += 1
    fasta_input_handle.close()

    return length_for_scaffold, scaffold_count

if __name__ == "__main__":
    main()