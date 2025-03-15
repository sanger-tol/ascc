#!/usr/bin/env python3

VERSION = "3.0.0"
DESCRIPTION = f"""
---
Script for creating a BlobToolKit dataset from the ASCC output files
Version: {VERSION}
---

Written by Eerik Aunin (ea10)

Modified by Damon-Lee Pointon (@dp24/@DLBPointon)

"""

import general_purpose_functions as gpf
import argparse
from pathlib import Path
import textwrap
import sys
import os.path
import re
from collections import defaultdict
from typing import Dict, List, Tuple


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="createBTKdatasets",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(DESCRIPTION),
    )
    parser.add_argument("-n", "--name", required=True, type=str, help="Assembly name (for the output files)")
    parser.add_argument(
        "-tn",
        "--taxon_name",
        required=True,
        type=str,
        help="The Taxon name of the assembly (Scientific name of the species + subspecies if applicable)",
    )
    parser.add_argument("-id", "--taxid", required=True, type=int, help="Taxon ID of the assembly")
    parser.add_argument(
        "-td", "--taxdump", required=True, type=str, help="Path to the directory containing the NCBI taxdump"
    )
    parser.add_argument("-f", "--fasta", required=True, type=str, help="The path for the assembly fasta file")
    parser.add_argument(
        "-d",
        "--dataset",
        type=str,
        required=True,
        help="The folder containing the data generated throughout the pipeline",
    )
    parser.add_argument("-bh", "--blastn_hits", default="N", type=str, help="Path to the BLASTN hits file")
    parser.add_argument(
        "-ud", "--uniprot_diamond_hits", default="N", type=str, help="Path to the UNIPROT diamond BlastX hits file"
    )
    parser.add_argument("-nr", "--nr_diamond_hits", default="N", type=str, help="Path to the DIAMOND BlastX hits file")
    parser.add_argument(
        "-r", "--mapped_reads", default="N", type=str, help="Path to mapped reads BAM for coverage estimation"
    )
    parser.add_argument("-t", "--tiara", default="N", type=str, help="Path to the tiara_out.txt file")
    parser.add_argument(
        "-p", "--pca", default="N", type=str, help="Path to the kmers_dim_reduction_embeddings_combined.csv file"
    )
    parser.add_argument("-fc", "--fcs_gx", default="N", type=str, help="Path to the fcs-gx_summary.csv.csv file")
    parser.add_argument("-k", "--kraken", default="N", type=str, help="Path to the nt_kraken_lineage.txt file")
    parser.add_argument("-ms", "--markerscan", default="N", type=str, help="Path to the cobiontid_markerscan.csv file")
    parser.add_argument("-cv", "--contigviz", default="N", type=str, help="Path to the contigviz_results.csv file")
    parser.add_argument("-o", "--output", default="btk_datasets", type=str, help="Output directory")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to utilise")
    parser.add_argument("--alias", type=str, default="", help="Assembly alias")
    parser.add_argument(
        "--dry_run", dest="dry_run", action="store_true", help="Dry run (print commands without executing)"
    )
    parser.add_argument("-v", "--version", action="version", version=VERSION)
    parser.add_argument("--debug", action="store_true", help="Print debug information")

    return parser.parse_args(argv)


def create_assembly_yaml(assembly_yaml_path, assembly_alias, taxon_name):
    """
    Creates the assembly YAML file for creating a BlobToolKit dataset
    """
    if ".gz" in assembly_alias:
        assembly_alias = assembly_alias.replace(".gz", "_gz")
    out_string = "assembly:\n  accession: NA\n  alias: {}\n  record_type: scaffold\n  bioproject: NA\n  biosample: NA\ntaxon:\n  name: {}".format(
        assembly_alias, taxon_name
    )
    with open(assembly_yaml_path, "w") as f:
        f.write(out_string)


def tiara_results_to_btk_format(tiara_results_path, outfile_path):
    """
    Reformatting Tiara output file so that the summarised results of the first and second pass of Tiara can be
        added to a BlobToolKit dataset
    """
    tiara_data = gpf.l(tiara_results_path)
    tiara_data = tiara_data[1 : len(tiara_data)]
    with open(outfile_path, "w") as f:
        f.write("identifier\ttiara\n")
        for line in tiara_data:
            split_line = line.split()
            if len(split_line) != 3:
                sys.stderr.write("Failed to parse the Tiara results file {}\n".format(tiara_results_path))
                sys.exit(1)
            first_pass_result = split_line[1]
            second_pass_result = split_line[2]
            if second_pass_result != "n/a":
                first_pass_result = second_pass_result
            f.write(split_line[0] + "\t" + first_pass_result + "\n")


def detect_dim_reduction_methods(kmers_dim_reduction_output_path, debug=False) -> Dict[str, Tuple[int, List[str], str]]:
    """
    Parses the header of the kmers dimensionality reduction report file to detect
    which dimensionality reduction methods were used and how many dimensions each has.
    Returns a dictionary where keys are sanitised method names and values are tuples of:
    (dimensions, column_names, original_method_name).
    """
    with open(kmers_dim_reduction_output_path) as f:
        header_string = f.readline().strip()

    split_header = header_string.split(",")

    # Get columns that start with embedding_dim_
    embedding_cols = [col for col in split_header if col.startswith("embedding_dim_")]

    if debug:
        sys.stderr.write(f"[DEBUG] Found embedding columns: {embedding_cols}\n")

    # Group columns by method
    method_columns = defaultdict(list)
    method_originals = {}  # Store original method names
    for col in embedding_cols:
        # Extract method name from the full column name
        method = re.sub(r"^embedding_dim_\d+_", "", col)
        # Create a sanitised version for the dictionary key
        sanitised_method = method.replace("-", "_")
        method_columns[sanitised_method].append(col)
        method_originals[sanitised_method] = method  # Store the original name

    # Sort columns for each method to ensure correct order
    for method in method_columns:
        method_columns[method].sort(key=lambda x: int(re.search(r"embedding_dim_(\d+)_", x).group(1)))

    # Create result dictionary with dimensions, sorted columns, and original method name
    result = {}
    for sanitised_method, columns in method_columns.items():
        dims = len(columns)
        original_method = method_originals[sanitised_method]
        result[sanitised_method] = (dims, columns, original_method)
        if debug:
            sys.stderr.write(
                f"[DEBUG] Detected method: {original_method} (sanitised: {sanitised_method}) with {dims} dimensions\n"
            )
            sys.stderr.write(f"[DEBUG] Columns: {columns}\n")

    return result


def shorten_method_name(method, debug=False):
    """
    Shortens specific method names for BlobToolKit variables
    """
    original = method
    if method.startswith("Autoencoder"):
        # Replace "Autoencoder" with "AE"
        method = method.replace("Autoencoder", "AE")
    elif method == "Non-Negative Matrix Factorization" or method == "Non_Negative_Matrix_Factorization":
        method = "NNMF"

    if debug and method != original:
        sys.stderr.write(f"[DEBUG] Shortened method name from '{original}' to '{method}'\n")
    return method


def sanitise_btk_variable(name, used_names=None, debug=False):
    """
    Sanitises variable names for BlobToolKit
    - Replaces hyphens with underscores
    - Truncates to 34 characters by default
    - If duplicates exist, truncates further to allow for suffix
    - Ensures total length never exceeds 36 characters
    """
    original = name
    name = name.replace("-", "_")

    # If no used_names provided or name not in used_names, truncate to 34 chars
    if not used_names or name[:34] not in [n[:34] for n in used_names]:
        name = name[:34]
        if debug and name != original:
            sys.stderr.write(f"[DEBUG] Sanitised name from '{original}' to '{name}' (no suffix needed)\n")
        return name

    # If we need to add a suffix, find the next available number
    base_name = name[:33]  # Leave room for _N
    suffix = 1
    while f"{base_name}_{suffix}" in used_names:
        suffix += 1
        # If suffix becomes double digits, shorten base_name further
        if suffix >= 10:
            base_name = name[:32]

    name = f"{base_name}_{suffix}"
    if debug:
        sys.stderr.write(f"[DEBUG] Sanitised name from '{original}' to '{name}' (with suffix)\n")
    return name


def set_default_plot_variables(args, command_list):
    """
    Sets the default plot variables for the BlobToolKit dataset based on available data
    """
    # Track available variables
    has_coverage = (
        args.mapped_reads != "N" and os.path.isfile(args.mapped_reads) and os.stat(args.mapped_reads).st_size > 0
    )
    has_kmers = args.pca != "N" and os.path.isfile(args.pca) and os.stat(args.pca).st_size > 0
    has_fcs_gx = args.fcs_gx != "N" and os.path.isfile(args.fcs_gx) and os.stat(args.fcs_gx).st_size > 0
    has_tiara = args.tiara != "N" and os.path.isfile(args.tiara) and os.stat(args.tiara).st_size > 0
    has_kraken = args.kraken != "N" and os.path.isfile(args.kraken) and os.stat(args.kraken).st_size > 0
    has_blast = any(
        n != "N" and os.path.isfile(n) and os.stat(n).st_size > 0
        for n in [args.blastn_hits, args.uniprot_diamond_hits, args.nr_diamond_hits]
    )

    # Determine X and Y variables
    x_var = "gc"  # Default X is always GC
    y_var = None

    if has_coverage:
        # Get the coverage variable name from the BAM filename
        bam_filename = os.path.basename(args.mapped_reads)
        coverage_var = os.path.splitext(bam_filename)[0] + "_cov"
        y_var = coverage_var
    elif has_kmers:
        # Use the first dimension of the first kmers method
        method_dimensions = detect_dim_reduction_methods(args.pca, args.debug)
        if method_dimensions:
            first_method = list(method_dimensions.keys())[0]
            shortened_method = shorten_method_name(first_method, args.debug)
            sanitised_method = sanitise_btk_variable(shortened_method, None, args.debug)
            x_var = f"embedding_dim_1_{sanitised_method}"
            y_var = f"embedding_dim_2_{sanitised_method}"
    elif has_fcs_gx:
        # Use fcs_gx_coverage_by_div as Y variable if no coverage or kmers are available
        y_var = "fcs_gx_coverage_by_div"

    # Determine categorical variable based on priority
    cat_var = None
    if has_fcs_gx:
        cat_var = "fcs_gx_div"
    elif has_tiara:
        cat_var = "tiara"
    elif has_kraken:
        cat_var = "nt_kraken_phylum"
    elif has_blast:
        cat_var = "bestsum_phylum"

    # Build the blobtools replace command
    replace_cmd = f"blobtools replace"
    if x_var:
        replace_cmd += f" --key plot.x={x_var}"
    if y_var:
        replace_cmd += f" --key plot.y={y_var}"
    if cat_var:
        replace_cmd += f" --key plot.cat={cat_var}"

    replace_cmd += f" {args.output}"

    # Add the command to the list if we have at least one key to set
    if "--key" in replace_cmd:
        command_list.append(replace_cmd)

    return command_list


def preprocess_fcsgx_csv(input_path, output_path, debug=False):
    """
    Preprocess the FCS-GX CSV file to replace "None" or empty values with "0" in numerical columns only.
    This helps BlobToolKit infer the correct data types for these columns.
    """
    if debug:
        sys.stderr.write(f"[DEBUG] Preprocessing FCS-GX CSV file: {input_path} -> {output_path}\n")

    # Read the input CSV file
    with open(input_path, "r") as f:
        lines = f.readlines()

    if not lines:
        sys.stderr.write(f"[ERROR] Empty or missing FCS-GX file: {input_path}\n")
        return

    # Parse the header to identify column indices
    header = lines[0].strip().split(",")

    # Identify indices of numerical columns
    numerical_columns = ["fcs_gx_coverage_by_div", "fcs_gx_coverage_by_tax", "fcs_gx_score"]
    numerical_indices = []

    for col in numerical_columns:
        if col in header:
            numerical_indices.append(header.index(col))
            if debug:
                sys.stderr.write(f"[DEBUG] Found numerical column '{col}' at index {header.index(col)}\n")
        else:
            sys.stderr.write(f"[WARNING] Column '{col}' not found in FCS-GX file header\n")

    # Process each line
    processed_lines = [lines[0]]  # Keep the header unchanged

    for line_num, line in enumerate(lines[1:], 1):  # Skip header
        fields = line.strip().split(",")

        # Ensure the line has enough fields
        if len(fields) < max(numerical_indices, default=-1) + 1:
            sys.stderr.write(f"[WARNING] Line {line_num} has fewer fields than expected, skipping\n")
            processed_lines.append(line)
            continue

        # Replace "None" or empty values with "0" only in numerical columns
        for idx in numerical_indices:
            if idx < len(fields):
                if fields[idx] == "None" or fields[idx] == "":
                    fields[idx] = "0"

        # Reconstruct the line
        processed_line = ",".join(fields) + "\n"
        processed_lines.append(processed_line)

    # Write the processed lines to the output file
    with open(output_path, "w") as f:
        f.writelines(processed_lines)

    if debug:
        sys.stderr.write(f"[DEBUG] FCS-GX CSV preprocessing complete\n")


def main(args):
    command_list = []

    assembly_alias = args.name if args.alias == "" else args.alias

    edited_assembly_title = args.name.replace(".", "_").replace(" ", "_")

    assembly_yaml_path = args.output + "/" + edited_assembly_title + "BTK_DS.yaml"

    if args.dry_run == False:
        Path(args.dataset).mkdir(parents=True, exist_ok=True)
        create_assembly_yaml(assembly_yaml_path, assembly_alias, args.taxon_name)

    # Base command for new BTK Dataset
    blobtools_create_command = f"blobtools create --fasta {args.fasta} --meta {assembly_yaml_path} --taxid {args.taxid} --taxdump {args.taxdump} {args.output}"
    gpf.run_system_command(blobtools_create_command, dry_run=args.dry_run)

    # ADDING BLAST HIT DATA TO BTK
    hits_file_paths = [args.blastn_hits, args.uniprot_diamond_hits, args.nr_diamond_hits]

    hits_file = [n for n in hits_file_paths if n != "N" and os.path.isfile(n) is True and os.stat(n).st_size > 0]

    if len(hits_file) > 0:
        add_hits_command = "blobtools add"
        for file in hits_file:  # Only include files that exist and are not empty
            add_hits_command += f" --hits {file}"
        add_hits_command += f" --taxrule bestsum --taxdump {args.taxdump} {args.output}"
        command_list.append(add_hits_command)

    # ADDING MAPPED READS DATA TO BTK
    if (
        args.mapped_reads != "N"
        and os.path.isfile(args.mapped_reads) is True
        and os.stat(args.mapped_reads).st_size > 0
    ):
        add_cov_command = f"blobtools add --cov {args.mapped_reads} --threads {args.threads} {args.output}"
        command_list.append(add_cov_command)

    # ADDING TIARA
    if args.tiara != "N" and os.path.isfile(args.tiara) and os.stat(args.tiara).st_size > 0:
        tiara_reformatted_output_path = args.dataset + "/tiara_out_btk_format.tsv"
        tiara_results_to_btk_format(args.tiara, tiara_reformatted_output_path)
        add_tiara_command = f"blobtools add --text {tiara_reformatted_output_path} --text-delimiter '\t' --text-cols 'identifier=identifiers,tiara=tiara' --text-header {args.output}"
        command_list.append(add_tiara_command)

    # ADDING KMER DIM REDUCTION
    if args.pca != "N" and os.path.isfile(args.pca) and os.stat(args.pca).st_size > 0:
        method_info = detect_dim_reduction_methods(args.pca, args.debug)
        used_names = set()
        for sanitised_method, (_, columns, _) in method_info.items():
            # Get the sanitised name for BlobToolKit variables
            shortened_method = shorten_method_name(sanitised_method, args.debug)
            sanitised_method = sanitise_btk_variable(shortened_method, used_names, args.debug)
            used_names.add(sanitised_method)

            # Create the text-cols string dynamically based on number of dimensions
            cols_list = ["scaff=identifiers"]
            # Use original column name from CSV (left side) and sanitised name for BTK variable (right side)
            for i, col in enumerate(columns, 1):
                btk_var = f"embedding_dim_{i}_{sanitised_method}"  # BTK variable name
                cols_list.append(f"{col}={btk_var}")  # Map original column to BTK variable
                if args.debug:
                    sys.stderr.write(f"[DEBUG] Column mapping: {col} -> {btk_var}\n")

            cols_string = ",".join(cols_list)

            add_embedding_command = (
                f"blobtools add --text {args.pca} "
                f"--text-delimiter ',' "
                f"--text-cols {cols_string} "
                f"--text-header {args.output}"
            )
            if args.debug:
                sys.stderr.write(f"[DEBUG] Generated blobtools command:\n{add_embedding_command}\n")
            command_list.append(add_embedding_command)

    # ADDING KRAKEN DATA
    if args.kraken != "N" and os.path.isfile(args.kraken) and os.stat(args.kraken).st_size > 0:
        for taxonomy_level in ("species", "genus", "family", "order", "class", "phylum", "kingdom", "domain"):
            add_kraken_command = f"blobtools add --text {args.kraken} --text-delimiter ',' --text-cols scaff=identifiers,nt_kraken_{taxonomy_level}=nt_kraken_{taxonomy_level} --text-header {args.output}"
            command_list.append(add_kraken_command)

    # ADDING FCS_GX DATA
    if args.fcs_gx != "N" and os.path.isfile(args.fcs_gx) and os.stat(args.fcs_gx).st_size > 0:
        # Preprocess the FCS-GX CSV file to ensure numerical columns are recognized as such
        fcsgx_preprocessed_path = args.dataset + "/fcs-gx_summary_preprocessed.csv"
        preprocess_fcsgx_csv(args.fcs_gx, fcsgx_preprocessed_path, args.debug)

        # Add categorical columns first
        add_fcs_gx_categorical_command = (
            f"blobtools add --text {fcsgx_preprocessed_path} "
            f"--text-delimiter ',' "
            f"--text-cols 'scaff=identifiers,fcs_gx_top_tax_name=fcs_gx_top_tax_name,fcs_gx_div=fcs_gx_div,fcs_gx_action=fcs_gx_action' "
            f"--text-header {args.output}"
        )
        if args.debug:
            sys.stderr.write(f"[DEBUG] Generated FCS-GX categorical command:\n{add_fcs_gx_categorical_command}\n")
        command_list.append(add_fcs_gx_categorical_command)

        # Add numerical columns separately
        add_fcs_gx_numerical_command = (
            f"blobtools add --text {fcsgx_preprocessed_path} "
            f"--text-delimiter ',' "
            f"--text-cols 'scaff=identifiers,fcs_gx_coverage_by_div=fcs_gx_coverage_by_div,"
            f"fcs_gx_coverage_by_tax=fcs_gx_coverage_by_tax,fcs_gx_score=fcs_gx_score' "
            f"--text-header {args.output}"
        )
        if args.debug:
            sys.stderr.write(f"[DEBUG] Generated FCS-GX numerical command:\n{add_fcs_gx_numerical_command}\n")
        command_list.append(add_fcs_gx_numerical_command)

    export_table_command = f"blobtools filter --table btk_summary_table_full.tsv {args.output}"
    command_list.append(export_table_command)

    # EXECUTE ALL BTK COMMANDS
    for i in command_list:
        gpf.run_system_command(i, dry_run=args.dry_run)

    # Set default plot variables
    default_var_commands = set_default_plot_variables(args, [])
    for cmd in default_var_commands:
        gpf.run_system_command(cmd, dry_run=args.dry_run)


if __name__ == "__main__":
    main(parse_args())
