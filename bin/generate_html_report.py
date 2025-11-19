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
        "--barcode_file", help="Published path to PacBio barcode results file"
    )
    parser.add_argument("--fcs_adaptor_euk_file", help="Published path to FCS-Adaptor eukaryotic results file")
    parser.add_argument("--fcs_adaptor_prok_file", help="Published path to FCS-Adaptor prokaryotic results file")
    parser.add_argument("--trim_ns_file", help="Published path to trim Ns results file")
    parser.add_argument(
        "--vecscreen_file", help="Published path to vecscreen results file"
    )
    parser.add_argument(
        "--autofilter_file", help="Published path to autofiltering results file"
    )
    parser.add_argument(
        "--merged_table_file", help="Published path to merged table results file"
    )
    parser.add_argument(
        "--phylum_coverage_file", help="Published path to phylum coverage data file"
    )
    parser.add_argument(
        "--kmers_dir", help="Staged path to k-mer profile results directory"
    )
    parser.add_argument("--fasta_sanitation_log", help="Published path to FASTA sanitation log file")
    parser.add_argument(
        "--fasta_length_filtering_log", help="Published path to FASTA length filtering log file"
    )
    parser.add_argument("--samplesheet", help="Published path to input samplesheet CSV file")
    parser.add_argument("--params_file", help="Published path to input parameters YAML file")
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
        "--outdir",
        help="Pipeline output directory (relative to launch_dir), e.g. 'test_01'",
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

    # Staged input paths to READ data from
    barcode_read = args.barcode_file if args.barcode_file else None
    fcs_euk_read = args.fcs_adaptor_euk_file if args.fcs_adaptor_euk_file else None
    fcs_prok_read = args.fcs_adaptor_prok_file if args.fcs_adaptor_prok_file else None
    trim_ns_read = args.trim_ns_file if args.trim_ns_file else None
    vecscreen_read = args.vecscreen_file if args.vecscreen_file else None
    autofilter_read = args.autofilter_file if args.autofilter_file else None
    merged_table_read = args.merged_table_file if args.merged_table_file else None
    fasta_sanitation_log_read = args.fasta_sanitation_log if args.fasta_sanitation_log else None
    fasta_length_filtering_log_read = args.fasta_length_filtering_log if args.fasta_length_filtering_log else None

    # Debug output
    print(f"Found barcode file (staged): {barcode_read}", file=sys.stderr)
    print(f"Found FCS euk file (staged): {fcs_euk_read}", file=sys.stderr)
    print(f"Found FCS prok file (staged): {fcs_prok_read}", file=sys.stderr)
    print(f"Found trim Ns file (staged): {trim_ns_read}", file=sys.stderr)
    print(f"Found vecscreen file (staged): {vecscreen_read}", file=sys.stderr)
    print(f"Found autofilter file (staged): {autofilter_read}", file=sys.stderr)
    print(f"Found merged table file (staged): {merged_table_read}", file=sys.stderr)
    print(
        f"Found FASTA sanitation log file: {fasta_sanitation_log_read}", file=sys.stderr
    )
    print(
        f"Found FASTA length filtering log file: {fasta_length_filtering_log_read}",
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

    # Get phylum coverage file path directly from argument
    phylum_coverage_file = args.phylum_coverage_file if args.phylum_coverage_file else None
    if phylum_coverage_file:
        print(f"Found phylum coverage file: {phylum_coverage_file}", file=sys.stderr)

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
    if fasta_sanitation_log_read or fasta_length_filtering_log_read:
        print(f"Processing FASTA sanitation log files", file=sys.stderr)
        fasta_sanitation_data = load_fasta_sanitation_log(
            fasta_sanitation_log_read, fasta_length_filtering_log_read
        )

    # Load data from files
    barcodes_check_data = (
        load_barcode_check_results(barcode_read)
        if barcode_read
        else "No barcode check results found."
    )
    fcs_adaptor_euk_report_data = (
        load_fcs_adaptor_results_as_table(fcs_euk_read)
        if fcs_euk_read
        else "No FCS-Adaptor eukaryotic results found."
    )
    fcs_adaptor_prok_report_data = (
        load_fcs_adaptor_results_as_table(fcs_prok_read)
        if fcs_prok_read
        else "No FCS-Adaptor prokaryotic results found."
    )
    trim_Ns_data = (
        load_trim_Ns_results(trim_ns_read)
        if trim_ns_read
        else "No trim Ns results found."
    )
    vecscreen_data = (
        load_vecscreen_results_as_table(vecscreen_read)
        if vecscreen_read
        else "No vecscreen results found."
    )
    autofiltering_data = (
        load_autofiltering_results_as_table(autofilter_read)
        if autofilter_read
        else "No autofiltering results found."
    )
    contamination_check_merged_table_data = (
        load_contamination_check_merged_table(merged_table_read)
        if merged_table_read
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
    phylum_coverage_read = args.phylum_coverage_file if args.phylum_coverage_file else None
    if phylum_coverage_read:
        print(
            f"Loading phylum coverage data from: {phylum_coverage_read}",
            file=sys.stderr,
        )
        phylum_coverage_data = load_phylum_coverage_data(phylum_coverage_read)

    # Helper to compute published display paths from staged inputs
    def publish_path(subdir, read_path):
        if not read_path:
            return None
        base = os.path.basename(read_path)
        if not base:
            return None
        if args.outdir and args.output_prefix:
            rel = os.path.join(args.outdir, args.output_prefix, subdir, base)
            return os.path.join(args.launch_dir, rel) if args.launch_dir else rel
        return None

    # Build display paths for templates (published locations)
    barcode_disp = publish_path("filter_barcode", barcode_read)
    fcs_euk_disp = publish_path("fcs_adaptor", fcs_euk_read)
    fcs_prok_disp = publish_path("fcs_adaptor", fcs_prok_read)
    trim_ns_disp = publish_path("trailingns", trim_ns_read)
    vecscreen_disp = publish_path("summarise_vecscreen_output", vecscreen_read)
    autofilter_disp = publish_path("autofilter", autofilter_read)
    merged_table_disp = publish_path("ascc_main_output", merged_table_read)
    phylum_coverage_disp = publish_path("ascc_main_output", phylum_coverage_read)
    reference_disp = publish_path("filtered_fasta", args.reference) if args.reference else None
    fasta_sanitation_disp = publish_path("filtered_fasta", fasta_sanitation_log_read)
    fasta_length_filtering_disp = publish_path("filtered_fasta", fasta_length_filtering_log_read)
    kmers_dir_disp = None
    if args.outdir and args.output_prefix:
        rel_dir = os.path.join(args.outdir, args.output_prefix, "kmer_data")
        kmers_dir_disp = os.path.join(args.launch_dir, rel_dir) if args.launch_dir else rel_dir

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
        # Source file display paths (published locations)
        reference_file=reference_disp,
        samplesheet_file=(params_dict.get("input") if params_dict and params_dict.get("input") else (args.samplesheet if args.samplesheet else None)),
        params_file=args.params_file if args.params_file else None,
        barcode_file=barcode_disp,
        fcs_adaptor_files=[p for p in [fcs_euk_disp, fcs_prok_disp] if p],
        trim_ns_file=trim_ns_disp,
        vecscreen_file=vecscreen_disp,
        autofilter_file=autofilter_disp,
        merged_table_file=merged_table_disp,
        phylum_coverage_file=phylum_coverage_disp,
        kmers_dir=kmers_dir_disp,
        fasta_sanitation_files=[p for p in [fasta_sanitation_disp, fasta_length_filtering_disp] if p],
        fcs_gx_report_file=publish_path("fcsgx_data", args.fcs_gx_report_txt) if args.fcs_gx_report_txt else None,
        fcs_gx_taxonomy_file=publish_path("fcsgx_data", args.fcs_gx_taxonomy_rpt) if args.fcs_gx_taxonomy_rpt else None,
    )

    # Render the HTML report and fail if not produced
    ok = render_html_report(data, args.template_dir, output_html)
    if not ok or not os.path.exists(output_html):
        print(f"Error: HTML report was not generated at {output_html}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
