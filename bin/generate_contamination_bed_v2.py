#!/usr/bin/env python
#
# -----------------------------
# Generate_contamination_bed
# -----------------------------
# This script creates a bed file summarising the comtamination found
# during the ASCC run. With the aim of informing the decontamination of
# the input fasta file.
# -----------------------------
# Adapted from ascc_nf_core_to_bed.py by James Torrance (jt8@sanger.ac.uk) by Will Eagles (we3@sanger.ac.uk)
# Re-written by Damon-Lee B Pointon (dp24@sanger.ac.uk)

import os
import sys
import re
import collections
import csv
import gzip
from Bio import SeqIO
import argparse

def parse_args(args=None):
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--assembly_name", type=str, help="Name of assembly for output", default="ASSEMBLY"
    )
    parser.add_argument(
        "--assembly_type", type=str, help="genomic or organellar", enum=["genomic", "organellar"]
    )
    parser.add_argument(
		"--assembly_path", type=str, help="The input assembly file (fasta)"
	)
    parser.add_argument(
        "--is_organelle", type=bool, help="Is the input assembly an organellar assembly?",
    )
    parser.add_argument(
        "--no_fcs_gx", type=bool, help="Do not use any FCS-GX results."
    )
    parser.add_argument(
        "--no_barcodes", type=bool, help="No multiplexing barcodes expected",
    )
    parser.add_argument(
        "--longread_paths", type=str, nargs="?", help="Paths of longreads used in ASCC run",
    )
    parser.add_argument("--version", action="version", version="2.0.0")

    return parser.parse_args(args)


def get_scaffold_for_fasta(assembly_file: str) -> tuple[dict[str, int], int]:
    """Get list of scaffolds from fasta file"""

    fasta_input_handle = gzip.open(assembly_file, "rt")
    length_for_scaffold = {}
    scaffold_count = 0
    for record in SeqIO.parse(fasta_input_handle, "fasta"):
        length_for_scaffold[record.id] = len(record.seq)
        scaffold_count += 1
    fasta_input_handle.close()

    return length_for_scaffold, scaffold_count


def parse_fcs_gx_file(fcs_gx_file: str) -> dict[str, str]:
    fcs_gx_taxonomy_for_scaffold = {}
    with open(fcs_gx_file) as fcs_gx_handle:
        next(fcs_gx_handle)  # Skip header
        for line in fcs_gx_handle:
            fields = re.split(",", line.rstrip())
            if len(fields) < 4:
                sys.exit(f"FCS-GX - Not enough fields in line {line}")
            scaffold = fields[0]
            taxonomy = fields[1]
            fcs_gx_div = fields[3]
            if re.search(":", fcs_gx_div):
                (fcs_gx_div_short, fcs_gx_div_long) = re.split(":", fcs_gx_div)
                taxonomy = f"{taxonomy}, {fcs_gx_div_long}"
            fcs_gx_taxonomy_for_scaffold[scaffold] = taxonomy

    return fcs_gx_taxonomy_for_scaffold


def parse_fcs_adaptor_file(input_file, index: dict[str, int]) -> dict[str, int]:
    sequences_removed = {}
    with open(input_file) as fcs_adaptor_handle:
        for line in fcs_adaptor_handle:
            fields = re.split(r"\s+", line.rstrip())
            if not re.match("^#", line):
                accession = fields[0]
                length = fields[1]
                action = fields[2]
                ranges = fields[3]
                start = 0
                end = 0
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
                    if accession in index:
                        end = index[accession]
                    else:
                        sys.exit(f"Cannot find scaffold {accession} in FASTA assembly")
                    treatment = "REMOVE"
                else:
                    sys.exit(f"Action not recognised: {action}")

                main_output = f"{accession}\t{str(start)}\t{str(end)}\t{treatment}\n"
                sequences_removed[main_output] = end - start
        return sequences_removed


def parse_barcode_files(barcode_files) -> dict[str, int]:
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


def parse_trim_ns_file(trim_ns_file) -> dict[str, int]:
    sequences_removed = {}
    with open(trim_ns_file) as trim_ns_handle:
        for line in trim_ns_handle:
            action, scaffold, start, end = re.split(r"\s+", line.rstrip())
            if action in (
                "TRIM",
                "REVCLIP",
                "FWDCLIP",
                "CLIP",
            ):
                treatment = "TRIM"
                main_output ="\t".join([scaffold, str(start), str(end), treatment]) + "\n"
                sequences_removed[main_output] = int(end) - int(start)
    return sequences_removed


def parse_organelle_csv(organelle_file: str):
    sequences_removed = {}
    scaffolds_removed = 0
    with open(organelle_file) as organelle_handle:
        bed_csv_reader = csv.reader( organelle_handle, delimiter="\t" )
        for field_set in bed_csv_reader:
            scaffold = field_set[0]
            start = int(field_set[1])
            end = int(field_set[2])
            (length, percentage) = re.split(",", field_set[3])
            percentage = float(percentage)
            # Discard anything over 95% coverage, or over 90% coverage if it's below 50 kb
            if ( percentage > 90 and end < 50000 ) or ( percentage > 95 ):
                treatment = "REMOVE"
                main_output = "\t".join([scaffold, str(start), str(end), treatment])+ "\n"
                sequences_removed[main_output] = end
    return sequences_removed, scaffolds_removed


def build_contamination_report(
    args.assembly_type,
    length_for_scaffold,
    lengths_removed,
    scaffolds_removed,
    scaffold_count,
    fcs_gx_taxonomy_removed,
    sequences_removed_by_reason,
    fcs_ambiguity_text
):
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


def main():
    # Collectors
    lengths_removed = []
    scaffolds_removed: int = 0
    fcs_gx_taxonomy_for_scaffold: dict[str, str] = {}
    sequences_removed_by_reason = collections.defaultdict(list)

    args = parse_args(None)

    files_for_decon = {
        "merged_filter_file": ( None if not args.no_fcs_gx and not args.is_organellar else args.abnormal_decon_csv ),
        "fcs_gx_file": ( None if not args.no_fcs_gx and not args.is_organellar else args.fcs_gx_file ),
        "fcs_adaptor_file": args.euk_fcs_adaptor_file,
        "trim_ns_file": (sys.exit("Trim Ns file is missing!") if args.trim_ns is None else args.trim_ns),
        "barcodes": (
            os.listdir(args.barcode_dir) if os.listdir(args.barcode_dir) > 1 else sys.exit("No barcode files")
        ),
        "organellar_recomendations": (
            os.listdir(args.organellar_recomendations_dir) if os.listdir(args.organellar_recomendations_dir) > 1 else []
        )
    }

    main_decon_bed_output_file: str = f"{args.assembly_name}.contamination.bed"
    tiara_decon_bed_output_file: str = f"{args.assembly_name}.tiara.bed"

    main_bed_output_handle = open(main_decon_bed_output_file, "w")
    tiara_bed_output_handle = open(tiara_decon_bed_output_file, "w")

    # Parse FASTA assembly
    length_for_scaffold, scaffold_count = get_scaffold_for_fasta(args.assembly_path)

    # Parse FCS-GX output
    fcs_gx_taxonomy_for_scaffold: dict[str, str] = ( {} if files_for_decon["fcs_gx_file"] == None else parse_fcs_gx_file(files_for_decon["fcs_gx_file"]) )

    # Parse FCS-Adaptor output
    fcs_sequences_removed: dict[str, int] = parse_fcs_adaptor_file(files_for_decon["fcs_adaptor_file"], length_for_scaffold)
    for key, val in fcs_sequences_removed.items():
        main_bed_output_handle.write(key)
        sequences_removed_by_reason["FCS-Adaptor"].append(val)

    # Parse Barcodes
    if not args.no_barcodes:
        main_bed_output_handle.write("# Barcodes" + "\n")
        sequences_for_removal = parse_barcode_files(files_for_decon["barcodes"])
        for key, val in sequences_for_removal.items():
            main_bed_output_handle.write(key)
            sequences_removed_by_reason["Barcodes"].append(val)
            lengths_removed.append(val)

    # Parse Trim Ns File
    main_bed_output_handle.write("# Trim Ns\n")
    trim_ns_sequences_removed = parse_trim_ns_file(files_for_decon["trim_ns_file"])
    for key, val in trim_ns_sequences_removed.items():
        main_bed_output_handle.write(key)
        sequences_removed_by_reason["Trim Ns"].append(val)
        lengths_removed.append(val)

    # Parse Organellar Recomendation Files
    long_organelle_name = {
        "mito": "Mitochondrion",
        "plastid": "Plastid",
    }

    if files_for_decon["organellar_recomendations"] != []:
        for organelle in ["mito", "plastid"]:
            for i in files_for_decon["organellar_recomendations"]:
                if organelle in i:
                    main_bed_output_handle.write(f"# {long_organelle_name[organelle]}\n")
                    organelle_sequences_removed, scaffolds_removed_for_organelle = parse_organelle_csv(organelle_file)
                    for key, val in organelle_sequences_removed.items():
                        main_bed_output_handle.write(key)
                        sequences_removed_by_reason[long_organelle_name[organelle]].append(val)
                        lengths_removed.append(val)
                    scaffolds_removed += scaffolds_removed_for_organelle

    # Process FCS-GX datsa
    fcs_gx_taxonomy_removed = collections.defaultdict(list)
    fcs_ambiguity_text = ""

    if files_for_decon["merged_filter_file"] != None:
        main_bed_output_handle.write("# Merged ASCC call\n")
        with open(files_for_decon["merged_filter_file"]) as merged_filter_handle:
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

if __name__ == "__main__":
    main()
