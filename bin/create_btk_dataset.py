#!/usr/bin/env python3

VERSION = "2.0.0"
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


def detect_dim_reduction_methods(kmers_dim_reduction_output_path):
    """
    Parses the header of the kmers dimensionality reduction report file to detect
    which dimensionality reduction methods were used and how many dimensions each has.
    Returns a dictionary where keys are method names and values are number of dimensions.

    The function extracts method names by removing the 'embedding_dim_X_' prefix from column names,
    preserving the complete method name including any underscores it may contain.
    """
    import re

    with open(kmers_dim_reduction_output_path) as f:
        header_string = f.readline().strip()

    split_header = header_string.split(",")
    dim_reduction_methods = {}

    # Get columns that start with embedding_dim_
    embedding_cols = [col for col in split_header if col.startswith("embedding_dim_")]

    # Extract method names by removing the embedding_dim_X_ prefix
    for col in embedding_cols:
        # Use regex to remove the prefix pattern 'embedding_dim_digits_'
        method = re.sub(r"^embedding_dim_\d+_", "", col)

        # Count dimensions for this method if we haven't seen it yet
        if method not in dim_reduction_methods:
            dims = sum(1 for c in embedding_cols if re.sub(r"^embedding_dim_\d+_", "", c) == method)
            dim_reduction_methods[method] = dims

    return dim_reduction_methods


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
        method_dimensions = detect_dim_reduction_methods(args.pca)
        for method, n_dims in method_dimensions.items():
            # Create the text-cols string dynamically based on number of dimensions
            cols_list = ["scaff=identifiers"]
            cols_list.extend([f"embedding_dim_{i}_{method}=embedding_dim_{i}_{method}" for i in range(1, n_dims + 1)])
            cols_string = ",".join(cols_list)

            add_embedding_command = (
                f"blobtools add --text {args.pca} "
                f"--text-delimiter ',' "
                f"--text-cols {cols_string} "
                f"--text-header {args.output}"
            )
            command_list.append(add_embedding_command)

    # ADDIND KRAKEN DATA
    if args.kraken != "N" and os.path.isfile(args.kraken) and os.stat(args.kraken).st_size > 0:
        for taxonomy_level in ("species", "genus", "family", "order", "class", "phylum", "kingdom", "domain"):
            add_kraken_command = f"blobtools add --text {args.kraken} --text-delimiter ',' --text-cols scaff=identifiers,nt_kraken_{taxonomy_level}=nt_kraken_{taxonomy_level} --text-header {args.output}"
            command_list.append(add_kraken_command)

    # ADDING FCS_GX DATA
    if args.fcs_gx != "N" and os.path.isfile(args.fcs_gx) and os.stat(args.fcs_gx).st_size > 0:
        add_fcs_gx_results_command = f"blobtools add --text {args.fcs_gx} --text-delimiter ',' --text-cols 'scaff=identifiers,fcs_gx_top_tax_name=fcs_gx_top_tax_name,fcs_gx_div=fcs_gx_div,fcs_gx_action=fcs_gx_action' --text-header {args.output}"
        command_list.append(add_fcs_gx_results_command)

    export_table_command = f"blobtools filter --table btk_summary_table_full.tsv {args.output}"
    command_list.append(export_table_command)

    # EXECUTE ALL BTK COMMANDS
    for i in command_list:
        gpf.run_system_command(i, dry_run=args.dry_run)


if __name__ == "__main__":
    main(parse_args())
