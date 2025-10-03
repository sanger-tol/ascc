#!/usr/bin/env python3

# Constants
DEFAULT_OUTPUT_PREFIX = "report"
SCRIPT_VERSION = "1.0"
DEFAULT_BTK_INCLUDED = "true"
NOT_PROVIDED_MSG = "Not provided"

DESCRIPTION = f"""
Script for generating the HTML report of the ASCC pipeline.
Version = {SCRIPT_VERSION}
Written by Eerik Aunin (@eeaunin)
"""

import argparse
import os
import sys
import pandas as pd


from html_report_loaders import (
    load_samplesheet,
    load_yaml_params,
    load_barcode_check_results,
    load_contamination_check_merged_table,
    load_phylum_coverage_data,
    load_fcs_adaptor_results_as_table,
    load_trim_Ns_results,
    load_vecscreen_results_as_table,
    load_autofiltering_results_as_table,
    load_fasta_sanitation_log,
    load_fcsgx_report_as_table,
    load_fcsgx_taxonomy_as_table,
    load_kmer_dim_reduction_results,
    process_reference_file_line_by_line,
    find_files_in_dir,
)
from html_report_renderer import render_html_report, prepare_report_data


def main():
    parser = argparse.ArgumentParser(
        description="Generate HTML report for ASCC pipeline results"
    )
    parser.add_argument(
        "--output_dir", required=True, help="Output directory for the report"
    )
    parser.add_argument(
        "--template_dir",
        required=True,
        help="Directory containing the Jinja2 templates",
    )
    parser.add_argument(
        "--reference", help="Reference genome file for assembly statistics"
    )
    parser.add_argument(
        "--barcode_dir", help="Directory containing PacBio barcode results"
    )
    parser.add_argument("--fcs_dir", help="Directory containing FCS-Adaptor results")
    parser.add_argument("--trim_ns_dir", help="Directory containing trim Ns results")
    parser.add_argument(
        "--vecscreen_dir", help="Directory containing vecscreen results"
    )
    parser.add_argument(
        "--autofilter_dir", help="Directory containing autofiltering results"
    )
    parser.add_argument(
        "--merged_dir", help="Directory containing merged table results"
    )
    parser.add_argument(
        "--coverage_dir", help="Directory containing phylum coverage data"
    )
    parser.add_argument(
        "--kmers_dir", help="Directory containing k-mer profile results"
    )
    parser.add_argument("--fasta_sanitation_log", help="FASTA sanitation log file")
    parser.add_argument(
        "--fasta_length_filtering_log", help="FASTA length filtering log file"
    )
    parser.add_argument("--samplesheet", help="Input samplesheet CSV file")
    parser.add_argument("--params_file", help="Input parameters YAML file")
    parser.add_argument(
        "--params_json", help="JSON string containing the params object"
    )
    parser.add_argument("--fcs_gx_report_txt", help="Path to FCS-GX report text file")
    parser.add_argument(
        "--fcs_gx_taxonomy_rpt", help="Path to FCS-GX taxonomy report file"
    )

    parser.add_argument(
        "--btk_published_path",
        help="Path to the published BlobToolKit dataset in the output directory",
    )
    parser.add_argument(
        "--btk_included",
        help="Whether BlobToolKit dataset creation was included in the run",
    )
    parser.add_argument(
        "--launch_dir", help="Directory from which the workflow was launched"
    )
    parser.add_argument(
        "--output_prefix",
        default=DEFAULT_OUTPUT_PREFIX,
        help="Prefix for the output HTML file",
    )
    parser.add_argument(
        "--pipeline_version", required=True, help="Version of the ASCC pipeline"
    )
    parser.add_argument("--version", action="version", version=SCRIPT_VERSION)

    args = parser.parse_args()

    # Set up output paths
    output_html = os.path.join(args.output_dir, "site", f"{args.output_prefix}.html")
    os.makedirs(os.path.dirname(output_html), exist_ok=True)

    # Find files in directories
    barcode_file = None
    if args.barcode_dir and os.path.exists(args.barcode_dir):
        barcode_files = find_files_in_dir(args.barcode_dir, extension=".txt")
        if barcode_files:
            barcode_file = barcode_files[0]

    fcs_euk_file = None
    fcs_prok_file = None
    if args.fcs_dir and os.path.exists(args.fcs_dir):
        fcs_files = find_files_in_dir(args.fcs_dir)
        for f in fcs_files:
            if "euk" in f.lower():
                fcs_euk_file = f
            elif "prok" in f.lower():
                fcs_prok_file = f

    trim_ns_file = None
    if args.trim_ns_dir and os.path.exists(args.trim_ns_dir):
        # Remove extension=".txt" to find any file
        trim_ns_files = find_files_in_dir(args.trim_ns_dir)
        if trim_ns_files:
            trim_ns_file = trim_ns_files[0]  # Take the first file found

    vecscreen_file = None
    if args.vecscreen_dir and os.path.exists(args.vecscreen_dir):
        # Remove extension=".txt" to find any file
        vecscreen_files = find_files_in_dir(args.vecscreen_dir)
        if vecscreen_files:
            vecscreen_file = vecscreen_files[0]  # Take the first file found

    autofilter_file = None
    if args.autofilter_dir and os.path.exists(args.autofilter_dir):
        autofilter_files = find_files_in_dir(args.autofilter_dir, extension=".csv")
        if autofilter_files:
            autofilter_file = autofilter_files[0]

    merged_table_file = None
    if args.merged_dir and os.path.exists(args.merged_dir):
        merged_files = find_files_in_dir(
            args.merged_dir, extension=".tsv"
        ) or find_files_in_dir(args.merged_dir, extension=".csv")
        if merged_files:
            merged_table_file = merged_files[0]

    # Find FASTA sanitation log files
    fasta_sanitation_log_file = None
    fasta_length_filtering_log_file = None
    if args.fasta_sanitation_log and os.path.exists(args.fasta_sanitation_log):
        fasta_sanitation_log_file = args.fasta_sanitation_log
    if args.fasta_length_filtering_log and os.path.exists(
        args.fasta_length_filtering_log
    ):
        fasta_length_filtering_log_file = args.fasta_length_filtering_log

    # Debug output
    print(f"Found barcode file: {barcode_file}", file=sys.stderr)
    print(f"Found FCS euk file: {fcs_euk_file}", file=sys.stderr)
    print(f"Found FCS prok file: {fcs_prok_file}", file=sys.stderr)
    print(f"Found trim Ns file: {trim_ns_file}", file=sys.stderr)
    print(f"Found vecscreen file: {vecscreen_file}", file=sys.stderr)
    print(f"Found autofilter file: {autofilter_file}", file=sys.stderr)
    print(f"Found merged table file: {merged_table_file}", file=sys.stderr)
    print(
        f"Found FASTA sanitation log file: {fasta_sanitation_log_file}", file=sys.stderr
    )
    print(
        f"Found FASTA length filtering log file: {fasta_length_filtering_log_file}",
        file=sys.stderr,
    )
    print(f"Kmers directory: {args.kmers_dir}", file=sys.stderr)

    # Process reference file if provided
    reference_summary = None
    if args.reference and os.path.exists(args.reference):
        print(f"Processing reference file: {args.reference}", file=sys.stderr)
        reference_summary = process_reference_file_line_by_line(args.reference)

    # Load samplesheet if provided
    samplesheet_data = None
    sample_name = None
    if args.samplesheet and os.path.exists(args.samplesheet):
        print(f"Processing samplesheet file: {args.samplesheet}", file=sys.stderr)
        samplesheet_data, sample_name = load_samplesheet(args.samplesheet)
        if sample_name:
            print(
                f"Extracted sample name from samplesheet: {sample_name}",
                file=sys.stderr,
            )
        else:
            print("Could not extract sample name from samplesheet", file=sys.stderr)

    # Look for phylum coverage data file after sample_name is set
    phylum_coverage_file = None
    # First check the coverage_dir if it exists
    if args.coverage_dir and os.path.exists(args.coverage_dir) and sample_name:
        phylum_coverage_pattern = f"{sample_name}_phylum_counts_and_coverage.csv"
        # First try to find a file specific to this sample
        phylum_coverage_files = [
            f
            for f in find_files_in_dir(args.coverage_dir, extension=".csv")
            if os.path.basename(f) == phylum_coverage_pattern
        ]

        # If not found, fall back to any phylum coverage file
        if not phylum_coverage_files:
            phylum_coverage_files = [
                f
                for f in find_files_in_dir(args.coverage_dir, extension=".csv")
                if os.path.basename(f).endswith("_phylum_counts_and_coverage.csv")
            ]

        if phylum_coverage_files:
            phylum_coverage_file = phylum_coverage_files[0]
            print(
                f"Found phylum coverage file in coverage_dir: {phylum_coverage_file}",
                file=sys.stderr,
            )

    # If not found in coverage_dir, fall back to merged_dir for backward compatibility
    if (
        not phylum_coverage_file
        and args.merged_dir
        and os.path.exists(args.merged_dir)
        and sample_name
    ):
        phylum_coverage_pattern = f"{sample_name}_phylum_counts_and_coverage.csv"
        # First try to find a file specific to this sample
        phylum_coverage_files = [
            f
            for f in find_files_in_dir(args.merged_dir, extension=".csv")
            if os.path.basename(f) == phylum_coverage_pattern
        ]

        # If not found, fall back to any phylum coverage file
        if not phylum_coverage_files:
            phylum_coverage_files = [
                f
                for f in find_files_in_dir(args.merged_dir, extension=".csv")
                if os.path.basename(f).endswith("_phylum_counts_and_coverage.csv")
            ]

        if phylum_coverage_files:
            phylum_coverage_file = phylum_coverage_files[0]
            print(
                f"Found phylum coverage file in merged_dir: {phylum_coverage_file}",
                file=sys.stderr,
            )

    # Load parameters from either JSON or YAML file
    yaml_params_data = None
    params_dict = {}
    kmer_length = None
    if args.params_json:
        print(f"Params JSON argument provided", file=sys.stderr)
        yaml_params_data, params_dict = load_yaml_params(params_json=args.params_json)
    elif args.params_file:
        print(f"YAML params file argument: '{args.params_file}'", file=sys.stderr)
        if os.path.exists(args.params_file):
            print(
                f"Processing YAML parameters file: {args.params_file}", file=sys.stderr
            )
            yaml_params_data, params_dict = load_yaml_params(file_path=args.params_file)
        else:
            print(
                f"YAML parameters file does not exist: {args.params_file}",
                file=sys.stderr,
            )
            # Try to find the file in the current directory
            print(f"Current working directory: {os.getcwd()}", file=sys.stderr)
            print(f"Directory contents: {os.listdir('.')}", file=sys.stderr)
    else:
        print("No parameters provided (neither JSON nor YAML file)", file=sys.stderr)

    # Extract k-mer length from parameters
    if params_dict and "kmer_length" in params_dict:
        kmer_length = params_dict["kmer_length"]
        print(f"Extracted k-mer length from parameters: {kmer_length}", file=sys.stderr)

    # Load FASTA sanitation logs if provided
    fasta_sanitation_data = None
    if fasta_sanitation_log_file or fasta_length_filtering_log_file:
        print(f"Processing FASTA sanitation log files", file=sys.stderr)
        fasta_sanitation_data = load_fasta_sanitation_log(
            fasta_sanitation_log_file, fasta_length_filtering_log_file
        )

    # Load data from files
    barcodes_check_data = (
        load_barcode_check_results(barcode_file)
        if barcode_file
        else "No barcode check results found."
    )
    fcs_adaptor_euk_report_data = (
        load_fcs_adaptor_results_as_table(fcs_euk_file)
        if fcs_euk_file
        else "No FCS-Adaptor eukaryotic results found."
    )
    fcs_adaptor_prok_report_data = (
        load_fcs_adaptor_results_as_table(fcs_prok_file)
        if fcs_prok_file
        else "No FCS-Adaptor prokaryotic results found."
    )
    trim_Ns_data = (
        load_trim_Ns_results(trim_ns_file)
        if trim_ns_file
        else "No trim Ns results found."
    )
    vecscreen_data = (
        load_vecscreen_results_as_table(vecscreen_file)
        if vecscreen_file
        else "No vecscreen results found."
    )
    autofiltering_data = (
        load_autofiltering_results_as_table(autofilter_file)
        if autofilter_file
        else "No autofiltering results found."
    )
    contamination_check_merged_table_data = (
        load_contamination_check_merged_table(merged_table_file)
        if merged_table_file
        else "No contamination check merged table found."
    )
    # Handle kmers directory more robustly
    kmers_results = None
    if args.kmers_dir:
        if os.path.exists(args.kmers_dir) and os.path.isdir(args.kmers_dir):
            # Check if directory is empty
            if os.listdir(args.kmers_dir):
                kmers_results = load_kmer_dim_reduction_results(args.kmers_dir)
            else:
                print(f"Kmers directory is empty: {args.kmers_dir}", file=sys.stderr)
        else:
            print(
                f"Kmers directory not found or not a directory: {args.kmers_dir}",
                file=sys.stderr,
            )

    # Parse FCS-GX report files into metadata and tables
    fcs_gx_report_metadata = None
    fcs_gx_report_table = None
    fcsgx_override_note = None
    if args.fcs_gx_report_txt and os.path.exists(args.fcs_gx_report_txt):
        print(
            f"Processing FCS-GX report file: {args.fcs_gx_report_txt}", file=sys.stderr
        )
        fcs_gx_report_metadata, fcs_gx_report_table = load_fcsgx_report_as_table(
            args.fcs_gx_report_txt
        )
        if fcs_gx_report_metadata:
            print(f"Successfully extracted FCS-GX report metadata", file=sys.stderr)
        if fcs_gx_report_table:
            print(
                f"Successfully converted FCS-GX report to HTML table", file=sys.stderr
            )

    fcs_gx_taxonomy_metadata = None
    fcs_gx_taxonomy_table = None
    if args.fcs_gx_taxonomy_rpt and os.path.exists(args.fcs_gx_taxonomy_rpt):
        print(
            f"Processing FCS-GX taxonomy file: {args.fcs_gx_taxonomy_rpt}",
            file=sys.stderr,
        )
        fcs_gx_taxonomy_metadata, fcs_gx_taxonomy_table = load_fcsgx_taxonomy_as_table(
            args.fcs_gx_taxonomy_rpt
        )
        if fcs_gx_taxonomy_metadata:
            print(f"Successfully extracted FCS-GX taxonomy metadata", file=sys.stderr)
        if fcs_gx_taxonomy_table:
            print(
                f"Successfully converted FCS-GX taxonomy to HTML table", file=sys.stderr
            )

    # Check for BTK Dataset
    btk_dataset_path = None
    btk_included = args.btk_included if args.btk_included else DEFAULT_BTK_INCLUDED

    # Store the published path and launch directory for debugging
    btk_published_path = (
        args.btk_published_path if args.btk_published_path else NOT_PROVIDED_MSG
    )
    launch_dir = args.launch_dir if args.launch_dir else NOT_PROVIDED_MSG
    print(f"Launch directory: {launch_dir}", file=sys.stderr)

    # Resolve the path using the launch directory if it's a relative path
    full_btk_path = None
    if args.btk_published_path:
        if not os.path.isabs(args.btk_published_path) and args.launch_dir:
            full_btk_path = os.path.join(args.launch_dir, args.btk_published_path)
            print(
                f"Resolving relative path '{args.btk_published_path}' using launch directory '{args.launch_dir}' to '{full_btk_path}'",
                file=sys.stderr,
            )
        else:
            full_btk_path = args.btk_published_path
            print(f"Using path as is: {full_btk_path}", file=sys.stderr)

        # Check if the resolved path exists
        if os.path.exists(full_btk_path) and os.path.isdir(full_btk_path):
            btk_dataset_path = full_btk_path
            print(f"Found BTK dataset at: {btk_dataset_path}", file=sys.stderr)
        else:
            btk_dataset_path = None
            print(f"BTK dataset not found at: {full_btk_path}", file=sys.stderr)
    else:
        btk_dataset_path = None
        print("BTK published path not provided.", file=sys.stderr)

    # Create meta object from output_prefix
    meta = {"id": args.output_prefix}

    # Load phylum coverage data if available
    phylum_coverage_data = None
    if phylum_coverage_file:
        print(
            f"Loading phylum coverage data from: {phylum_coverage_file}",
            file=sys.stderr,
        )
        phylum_coverage_data = load_phylum_coverage_data(phylum_coverage_file)

    # Determine fcs_override flag from params if available and set note when raw files are missing
    fcs_override_flag = False
    try:
        if params_dict and isinstance(params_dict, dict):
            fcs_override_flag = bool(params_dict.get("fcs_override", False))
    except Exception:
        fcs_override_flag = False

    if fcs_override_flag and not (fcs_gx_report_metadata or fcs_gx_report_table or fcs_gx_taxonomy_metadata):
        fcsgx_override_note = (
            "FCS-GX was run externally (--fcs_override). Full raw report files were not provided to the pipeline, so the detailed FCS-GX tabs are omitted."
        )

    # Prepare data for the report
    data = prepare_report_data(
        reference_summary=reference_summary,
        samplesheet_data=samplesheet_data,
        sample_name=sample_name,
        yaml_params_data=yaml_params_data,
        barcodes_check_data=barcodes_check_data,
        fcs_adaptor_euk_report_data=fcs_adaptor_euk_report_data,
        fcs_adaptor_prok_report_data=fcs_adaptor_prok_report_data,
        trim_Ns_data=trim_Ns_data,
        vecscreen_data=vecscreen_data,
        autofiltering_data=autofiltering_data,
        contamination_check_merged_table_data=contamination_check_merged_table_data,
        coverage_per_phylum_data=phylum_coverage_data,
        kmers_results=kmers_results,
        fasta_sanitation_data=fasta_sanitation_data,
        # FCS-GX data
        fcs_gx_report_metadata=fcs_gx_report_metadata,
        fcs_gx_report_table=fcs_gx_report_table,
        fcs_gx_taxonomy_metadata=fcs_gx_taxonomy_metadata,
        fcs_gx_taxonomy_table=fcs_gx_taxonomy_table,
        fcsgx_override_note=fcsgx_override_note,
        timestamp=pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
        version=args.pipeline_version,  # Use the passed pipeline version
        meta=meta,  # Pass meta object
        kmer_length=kmer_length,  # Pass k-mer length
        btk_dataset_path=btk_dataset_path,  # Pass BlobToolKit dataset path
        btk_included=btk_included,  # Pass BlobToolKit included flag
        btk_published_path=btk_published_path,  # Pass BlobToolKit published path for debugging
        launch_dir=launch_dir,  # Pass launch directory for debugging
    )

    # Render the HTML report
    render_html_report(data, args.template_dir, output_html)


if __name__ == "__main__":
    main()
