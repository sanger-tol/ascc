#!/usr/bin/env python3

import os
import sys
from jinja2 import Environment, FileSystemLoader


def render_html_report(data, template_dir, output_file):
    """Render the HTML report using Jinja2 templates.

    Args:
        data (dict): Data to pass to the template
        template_dir (str): Directory containing the Jinja2 templates
        output_file (str): Path to the output HTML file
    """
    try:
        # Set up Jinja2 environment
        env = Environment(
            loader=FileSystemLoader(
                [template_dir, os.path.join(os.environ.get("WORKFLOW_PROJECTDIR", "."), "assets/templates")]
            )
        )

        # Get the template
        template = env.get_template("ascc_report.html.jinja")

        # Render template with data
        html_content = template.render(**data)

        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

        # Write HTML to output file
        with open(output_file, "w") as f:
            f.write(html_content)

        print(f"HTML report generated: {output_file}")
        return True
    except Exception as e:
        print(f"Error rendering HTML report: {e}", file=sys.stderr)
        return False


def prepare_report_data(
    reference_summary=None,
    samplesheet_data=None,
    sample_name=None,
    yaml_params_data=None,
    barcodes_check_data="No barcode check results found.",
    fcs_adaptor_euk_report_data="No FCS-Adaptor eukaryotic results found.",
    fcs_adaptor_prok_report_data="No FCS-Adaptor prokaryotic results found.",
    trim_Ns_data="No trim Ns results found.",
    vecscreen_data="No vecscreen results found.",
    autofiltering_data="No autofiltering results found.",
    contamination_check_merged_table_data="No contamination check merged table found.",
    coverage_per_phylum_data=None,
    kmers_results=None,
    fasta_sanitation_data=None,
    # FCS-GX report data
    fcs_gx_report_content=None,         # Original FCS-GX report content arg
    fcs_gx_taxonomy_content=None,       # Original FCS-GX taxonomy content arg
    fcs_gx_report_metadata=None,        # New FCS-GX report metadata arg
    fcs_gx_report_table=None,           # New FCS-GX report table arg
    fcs_gx_taxonomy_metadata=None,      # New FCS-GX taxonomy metadata arg
    fcs_gx_taxonomy_table=None,         # New FCS-GX taxonomy table arg
    timestamp=None,
    version="1.0",
    meta=None,                          # Meta parameter
    kmer_length=None,                   # K-mer length parameter
    btk_dataset_path=None,              # BlobToolKit dataset path parameter
    btk_included=None,                  # Whether BlobToolKit dataset creation was included in the run
    btk_published_path=None,            # Path being checked for the BlobToolKit dataset
    launch_dir=None,                    # Directory from which the workflow was launched
):
    """Prepare data for the HTML report.

    Args:
        reference_summary: HTML table with assembly statistics
        samplesheet_data: HTML table with samplesheet data
        sample_name: Name of the sample
        yaml_params_data: HTML table with YAML parameters
        barcodes_check_data: PacBio barcode check results
        fcs_adaptor_euk_report_data: FCS adaptor eukaryotic results
        fcs_adaptor_prok_report_data: FCS adaptor prokaryotic results
        trim_Ns_data: Trim Ns results
        vecscreen_data: VecScreen results
        autofiltering_data: Autofiltering results
        contamination_check_merged_table_data: Contamination check merged table
        coverage_per_phylum_data: Coverage per phylum data
        kmers_results: K-mer dimensionality reduction results
        fasta_sanitation_data: FASTA sanitation log data
        fcs_gx_report_content: Content of the FCS-GX report text file
        fcs_gx_taxonomy_content: Content of the FCS-GX taxonomy report file
        timestamp: Timestamp for the report
        version: Version of the pipeline

    Returns:
        dict: Data for the HTML report
    """
    # Use current timestamp if not provided
    if timestamp is None:
        import pandas as pd

        timestamp = pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")

    # Prepare data dictionary
    data = {
        "reference_summary": reference_summary,
        "samplesheet_data": samplesheet_data,
        "sample_name": sample_name,
        "yaml_params_data": yaml_params_data,
        "barcodes_check_data": barcodes_check_data,
        "fcs_adaptor_euk_report_data": fcs_adaptor_euk_report_data,
        "fcs_adaptor_prok_report_data": fcs_adaptor_prok_report_data,
        "trim_Ns_data": trim_Ns_data,
        "vecscreen_data": vecscreen_data,
        "autofiltering_data": autofiltering_data,
        "contamination_check_merged_table_data": contamination_check_merged_table_data,
        "coverage_per_phylum_data": coverage_per_phylum_data,
        "kmers_results": kmers_results,
        "fasta_sanitation_data": fasta_sanitation_data,
        # Original FCS-GX content (for backward compatibility)
        "fcs_gx_report_content": fcs_gx_report_content,
        "fcs_gx_taxonomy_content": fcs_gx_taxonomy_content,
        # New formatted FCS-GX data
        "fcs_gx_report_metadata": fcs_gx_report_metadata,
        "fcs_gx_report_table": fcs_gx_report_table,
        "fcs_gx_taxonomy_metadata": fcs_gx_taxonomy_metadata,
        "fcs_gx_taxonomy_table": fcs_gx_taxonomy_table,
        "timestamp": timestamp,
        "version": version,
        "meta": meta,
        "kmer_length": kmer_length,  # Add k-mer length to the data dictionary
        "btk_dataset_path": btk_dataset_path,  # Add BlobToolKit dataset path to the data dictionary
        "btk_included": btk_included,  # Add BlobToolKit included flag to the data dictionary
        "btk_published_path": btk_published_path,  # Add BlobToolKit published path for debugging
        "launch_dir": launch_dir,  # Add launch directory for debugging
    }

    return data
