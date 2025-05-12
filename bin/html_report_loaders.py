#!/usr/bin/env python3
"""
HTML report loaders for ASCC.

This module imports all the report loader functions from the report_loaders package.
The functions are used to load and format various types of data for inclusion in the HTML report.

Note: This file is a wrapper around the report_loaders package to maintain backward compatibility.
New code should import directly from the report_loaders package.
"""

# Import all loaders from the report_loaders package
from report_loaders import (
    # Basic loaders
    load_samplesheet,
    load_yaml_params,
    load_barcode_check_results,
    
    # Table loaders
    load_contamination_check_merged_table,
    load_phylum_coverage_data,
    load_trim_Ns_results,
    
    # Text loaders
    load_fcs_adaptor_results,
    load_vecscreen_results,
    load_autofiltering_results,
    
    # FASTA loaders
    load_fasta_length_filtering_log,
    load_fasta_sanitation_log,
    
    # FCS-GX loaders
    parse_fcsgx_metadata,
    flatten_nested_dict,
    metadata_to_html_table,
    format_fcsgx_metadata,
    fallback_load_fcsgx_taxonomy,
    load_fcsgx_taxonomy_as_table,
    
    # FCS-GX report loaders
    load_fcsgx_results,
    load_fcsgx_report_as_table,
    
    # K-mer loaders
    load_kmer_dim_reduction_results,
    
    # Reference loaders
    process_reference_file_line_by_line,
    
    # Utility functions
    find_files_in_dir,
)

# Make all loaders available at the module level
__all__ = [
    # Basic loaders
    'load_samplesheet',
    'load_yaml_params',
    'load_barcode_check_results',
    
    # Table loaders
    'load_contamination_check_merged_table',
    'load_phylum_coverage_data',
    'load_trim_Ns_results',
    
    # Text loaders
    'load_fcs_adaptor_results',
    'load_vecscreen_results',
    'load_autofiltering_results',
    
    # FASTA loaders
    'load_fasta_length_filtering_log',
    'load_fasta_sanitation_log',
    
    # FCS-GX loaders
    'parse_fcsgx_metadata',
    'flatten_nested_dict',
    'metadata_to_html_table',
    'format_fcsgx_metadata',
    'fallback_load_fcsgx_taxonomy',
    'load_fcsgx_taxonomy_as_table',
    
    # FCS-GX report loaders
    'load_fcsgx_results',
    'load_fcsgx_report_as_table',
    
    # K-mer loaders
    'load_kmer_dim_reduction_results',
    
    # Reference loaders
    'process_reference_file_line_by_line',
    
    # Utility functions
    'find_files_in_dir',
]
