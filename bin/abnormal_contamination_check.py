#!/usr/bin/env python

import general_purpose_functions as gpf
import sys
import os.path
import pathlib
import argparse
import textwrap

VERSION = "V1.1.0"

DESCRIPTION = """
-------------------------------------
    Abnormal Contamination Check
        Version = {VERSION}
-------------------------------------
Written by James Torrance
Modified by Eerik Aunin
Modified by Damon-Lee Pointon
-------------------------------------

Script for determining if there is
enough contamination found by FCS-GX
to warrant an abnormal contamination
report alarm. Partially based on code
written by James Torrance
-------------------------------------

"""


def parse_args():
    parser = argparse.ArgumentParser(
        prog="Abnormal Contamination Check",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(DESCRIPTION),
    )
    parser.add_argument("assembly", type=str, help="Path to the fasta assembly file")
    parser.add_argument("summary_path", type=str, help="Path to the tiara summary file")
    parser.add_argument("-q", "--out_prefix", type=str, help="Output file prefix for the report")
    parser.add_argument("-o", "--output", type=str, help="Path to output file", default="fcs-gx_alarm_indicator_file.txt")
    parser.add_argument("-p", "--alarm_percentage", type=int, help="Percentage of putative contaminant sequence in genomic assembly that will trip the alarm", default=3)
    parser.add_argument("-l", "--alarm_length_removed", type=int, help="Length of removed sequence is greater than default, greater than this will trip the alarm.", default=1e7)
    parser.add_argument("-s", "--alarm_scaff_length", type=int, help="Length of largest scaffold removed to trip alarm.", default=1.8e6)
    parser.add_argument("-t", "--alarm_scaff_percent_removed", type=float, help="Percentage of Scaffolds set for removal from assembly to trip the alarm.", default=10.0)
    parser.add_argument("-v", "--version", action="version", version=VERSION)
    return parser.parse_args()


def get_sequence_lengths(assembly_fasta_path):
    """
    Gets sequence lengths of a FASTA file and returns them as a dictionary
    """
    seq_lengths_dict = dict()
    fasta_data = gpf.read_fasta_in_chunks(assembly_fasta_path)
    for header, seq in fasta_data:
        seq_len = len(seq)
        seq_lengths_dict[header] = dict()
        seq_lengths_dict[header]["seq_len"] = seq_len
    return seq_lengths_dict


def load_fcs_gx_results(seq_dict, fcs_gx_and_tiara_summary_path):
    """
    Loads FCS-GX actions from the FCS-GX and Tiara results summary file, adds them to the dictionary that contains sequence lengths
    """
    fcs_gx_and_tiara_summary_data = gpf.l(fcs_gx_and_tiara_summary_path)
    fcs_gx_and_tiara_summary_data = fcs_gx_and_tiara_summary_data[
        1 : len(fcs_gx_and_tiara_summary_data)
    ]
    for line in fcs_gx_and_tiara_summary_data:
        split_line = line.split(",")
        assert len(split_line) == 5
        seq_name = split_line[0]
        fcs_gx_action = split_line[1]
        seq_dict[seq_name]["fcs_gx_action"] = fcs_gx_action
    return seq_dict


def main():
    args = parse_args()
    if os.path.isfile(args.summary_path) is False:
        sys.stderr.write(
            f"The FCS-GX and Tiara results file was not found at the expected location ({args.summary_path})\n"
        )
        sys.exit(1)

    if os.path.isfile(args.assembly) is False:
        sys.stderr.write(
            f"The assembly FASTA file was not found at the expected location ({args.assembly})\n"
        )
        sys.exit(1)

    seq_dict = get_sequence_lengths(args.assembly)
    seq_dict = load_fcs_gx_results(seq_dict, args.summary_path)

    total_assembly_length = 0
    lengths_removed = list()
    scaffolds_removed = 0
    scaffold_count = len(seq_dict)

    for seq_name in seq_dict:
        seq_len = seq_dict[seq_name]["seq_len"]
        if seq_dict[seq_name]["fcs_gx_action"] == "EXCLUDE":
            lengths_removed.append(seq_len)
            scaffolds_removed += 1
        total_assembly_length += seq_len

    alarm_threshold_for_parameter = {
        "TOTAL_LENGTH_REMOVED": args.alarm_length_removed,
        "PERCENTAGE_LENGTH_REMOVED": args.alarm_percentage,
        "LARGEST_SCAFFOLD_REMOVED": args.alarm_scaff_length,
        "PERCENTAGE_SCAFFOLDS_REMOVED": args.alarm_scaff_percent_removed
    }

    report_dict = {
        "TOTAL_LENGTH_REMOVED": sum(lengths_removed),
        "PERCENTAGE_LENGTH_REMOVED": 100 * sum(lengths_removed) / total_assembly_length,
        "LARGEST_SCAFFOLD_REMOVED": max(lengths_removed, default=0),
        "SCAFFOLDS_REMOVED": scaffolds_removed,
        "PERCENTAGE_SCAFFOLDS_REMOVED": 100 * scaffolds_removed / scaffold_count,
    }

    # Seperated out to ensure that the file is written in one go and doesn't confuse Nextflow
    with open(f"{args.out_prefix}_raw_report.txt", "a") as f:
        for x, y in report_dict.items():
            f.write(": ".join([x, str(y)]) + "\n")

    pathlib.Path(args.output).unlink(missing_ok=True)

    alarm_list = []
    stage1_decon_pass_flag = True
    for param in alarm_threshold_for_parameter:
        param_value = report_dict[param]
        alarm_threshold = alarm_threshold_for_parameter[param]

        # IF CONTAMINATING SEQ FOUND FILL FILE WITH ABNORMAL CONTAM
        if param_value > alarm_threshold_for_parameter[param]:
            stage1_decon_pass_flag = False
            alarm_list.append(
                f"YES_ABNORMAL_CONTAMINATION: Stage 1 decon alarm triggered for {param}: the value for this parameter in this assembly is {param_value} | alarm threshold is {alarm_threshold}\n"
            )

    # Seperated out to ensure that the file is written in one go and doesn't confuse Nextflow
    with open(args.output, "a") as f:
        f.write("".join(alarm_list))

    # IF NO CONTAM FILL FILE WITH NO CONTAM
    if stage1_decon_pass_flag is True:
        alarm_message = f"NO_ABNORMAL_CONTAMINATION: No scaffolds were tagged for removal by FCS-GX\n"
        with open(args.output, "a") as f:
            f.write(alarm_message)


if __name__ == "__main__":
    main()
