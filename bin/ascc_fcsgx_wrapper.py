#!/usr/bin/env python3

VERSION = "1.0.0"
DESCRIPTION = f"""
ASCC-FCSGX-WRAPPER.py
Version: {VERSION}
----------------------
By Damon-Lee Pointon (dp24)

Wrapper has been written to get arround an
issue with the singularity container of FCS not
working well on LSF.

This script runs the:
   - FILTER_FASTA_BY_LENGTH
    - Not running this will result in less contamination
        due to how fcs uses the assembly.

   - FCS-GX binary
    - The fcs program

   - PARSE_FCS_OUTPUT
    - Takes the fcs output directory and parses it into a
        singular file for use in ASCC

This script will then generate a fcs_samplesheet.csv containing:

    assembly,assembly_type,fcs_file
    iyTipFemo1,PRIMARY,*fcs_summary.csv
"""

FILTER_FA_SCRIPT="/lustre/scratch124/tol/teams/tola/users/dp24/ascc/bin/filter_fasta_by_length.py"
RUN_FCSGX_BINARY="/lustre/scratch124/tol/teams/tola/users/dp24/fcs-gx/dist/run_gx"
RUN_FCSGX_DATABS="/tmp/tol_data/fcs-gx/2023-01-24/"
PARSE_FCS_SCRIPT="/lustre/scratch124/tol/teams/tola/users/dp24/ascc/bin/parse_fcsgx_result.py"

# Eeriks module for general functions
import general_purpose_functions as gpf

import subprocess
import argparse
import textwrap
import sys
import os



def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="ASCC-FCS-WRAPPER",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(DESCRIPTION),
    )
    parser.add_argument(
        "-i",
        "--input_ascc_samplesheet",
        type=str,
        help="Path to ASCC input samplesheet"
    )
    parser.add_argument(
        "-n",
        "--ncbi_taxonomy_path",
        type=str,
        help="Path to NCBI taxonomy DB",
        default="/data/tol/resources/taxonomy/latest/new_taxdump/rankedlineage.dmp"
    )
    parser.add_argument(
        "-t",
        "--taxid",
        type=str,
        help="Taxid for assembly"
    )

    # FILTER_FASTA SPECIFIC ARGS
    parser.add_argument(
        "-c",
        "--cutoff",
        type=int,
        default=100000000,
        help="Cutoff value for filtering"
    )
    parser.add_argument(
        "-l",
        "--low_pass",
        action="store_true",
        help="Optional: low pass filtering mode (sequences longer than the cutoff value will be removed)",
    )
    parser.add_argument(
        "--remove_original_fasta",
        action="store_true",
        help="Optional: remove the input FASTA file after creating the filtered FASTA file",
    )

    parser.add_argument("-v", action="version", version=VERSION)

    return parser.parse_args()


def read_samplesheet(input_samplesheet: str) -> dict:
    output_dict = dict()

    print(f"PARSING ASCC SAMPLESHEET: {input_samplesheet}")
    with open(input_samplesheet) as input:
        for line in input:
            if line.startswith("sample"):
                pass
            else:
                split_line = line.split(",")
                assembly_name_full = f"{split_line[0]}_{split_line[1]}"
                fasta = split_line[2]
                output_dict[assembly_name_full] = fasta.strip()

    return output_dict


def gunzip_input(samplesheet: dict) -> dict:
    new_dict = dict()

    for x, y in samplesheet.items():
        if os.path.exists(f"./{x}.fasta"):
            print(f"{x}: ALREADY LOCAL AND UNZIPPED")
            new_file_path = os.path.abspath(f"./{x}.fasta")
            new_dict[x] = new_file_path
        else:
            if y.endswith(".gz"):
                print(f"{x}: COPYING TO ./")
                copy_args = ["cp", y, f"./{x}.fasta.gz"]
                subprocess.call(copy_args)

                print(f"{x}: UNZIPPING")
                gunzip_args = ["gunzip", f"./{x}.fasta.gz"]
                subprocess.call(gunzip_args)

                new_file_path = os.path.abspath(f"./{x}.fasta.gz")
                new_dict[x] = new_file_path
            else:
                sys.exit(1, "A file probably ends with fa rather than fasta(.gz)")

    return new_dict


def filter_fasta_by_length_wrapper(input_dict: dict, cutoff: int) -> dict:
    filter_dict = dict()
    for x, y in input_dict.items():
        if os.path.exists(f"{x}_filtering/{x}_filtered.fasta"):
            print(f"{x}: ALREADY FILTERED")
            filter_dict[x] = str(os.path.abspath(f"{x}_filtering/{x}_filtered.fasta"))
        else:
            print(f"{x}: FILTERING")
            mkdir_args = ["mkdir", "-p", f"{x}_filtering"]
            subprocess.call(mkdir_args)

            filter_args = [
                "python3",f"{FILTER_FA_SCRIPT}",
                str(y), str(cutoff), "--low_pass",
                "--remove_original_fasta",
                ">", f"{x}_filtering/{x}_filtered.fasta"
            ]
            subprocess.run(" ".join(filter_args), shell=True)
            filter_dict[x] = str(os.path.abspath(f"{x}_filtering/{x}_filtered.fasta"))

    return filter_dict


def run_fcs_gx(input_dict: dict, taxid) -> dict:
    fcs_dict = dict()
    for x, y in input_dict.items():
        print(f"{x}: RUN FCS_GX : STARTING")
        mkdir_args = ["mkdir", "-p", f"{x}_fcs"]
        subprocess.call(mkdir_args)

        fcs_args = [
            RUN_FCSGX_BINARY,
            "--fasta", y,
            "--gx-db", RUN_FCSGX_DATABS,
            "--tax-id", taxid,
            "--generate-logfile", "true",
            "--out-basename", x,
            "--out-dir", f"{x}_fcs"
        ]
        subprocess.run(" ".join(fcs_args), shell=True)

        # we need the fcs output folder
        fcs_dict[x] = str(os.path.abspath(f"./{x}_fcs"))
        print(f"{x}: RUN FCS_GX : COMPLETE")

    return fcs_dict


def parse_fcs_results(input_dict, ncbi_tax_path) -> dict:
    parsed_dict = dict()

    for x, y in input_dict.items():
        print(f"{x}: PARSING")
        mkdir_args = ["mkdir", "-p", f"{x}_filtering"]
        subprocess.call(mkdir_args)

        filter_args = [
            "python3",f"{PARSE_FCS_SCRIPT}",
            y, ncbi_tax_path,
            ">", f"{x}_parsed_fcsgx.csv"
        ]
        subprocess.run(" ".join(filter_args), shell=True)
        parsed_dict[x] = str(os.path.abspath(f"{x}_parsed_fcsgx.csv"))

    return parsed_dict


def generate_fcs_samplesheet(input_dict: dict) -> None:

    print("WRITING SAMPLESHEET")
    with open("fcs_samplesheet.csv", "w") as output:
        output.write("assembly,assembly_type,fcs_file")
        for x, y in input_dict:
            split_name = x.split("_")
            output.write(f"{split_name[0]},{split_name[1]},y")
    samplesheet = str(os.path.abspath("fcs_samplesheet.csv"))

    print(f"FINISHED : SAMPLESHEET -> {samplesheet}")


def main():
    args = parse_args()

    sample_dictionary = read_samplesheet(args.input_ascc_samplesheet)

    gunzipped_dictionary = gunzip_input(sample_dictionary)

    filtered_dict = filter_fasta_by_length_wrapper(gunzipped_dictionary, args.cutoff)
    print(filtered_dict)

    fcs_dict = run_fcs_gx(filtered_dict, args.taxid)
    print(fcs_dict)

    parsed_fcs_dict = parse_fcs_results(fcs_dict, args.ncbi_taxonomy_path)
    print(parsed_fcs_dict)

    generate_fcs_samplesheet(parsed_fcs_dict)

if __name__ == "__main__":
    main()
