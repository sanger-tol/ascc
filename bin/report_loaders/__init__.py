"""
Report loaders package for ASCC HTML report generation.

This package contains modules for loading and formatting various types of data
for inclusion in the HTML report.
"""

# Import commonly used functions for easier access
from .fcsgx_loaders import format_fcsgx_metadata
from .fcsgx_taxonomy_fallback import fallback_load_fcsgx_taxonomy

# Basic loaders
from .basic_loaders import (
    load_samplesheet,
    load_yaml_params,
    load_barcode_check_results,
)

# Table loaders
from .table_loaders import (
    load_contamination_check_merged_table,
    load_trim_Ns_results,
)

# Text loaders
from .text_loaders import (
    load_fcs_adaptor_results,
    load_vecscreen_results,
    load_autofiltering_results,
)

# FASTA loaders
from .fasta_loaders import (
    load_fasta_length_filtering_log,
    load_fasta_sanitation_log,
)

# FCS-GX loaders
from .fcsgx_loaders import (
    parse_fcsgx_metadata,
    flatten_nested_dict,
    metadata_to_html_table,
    format_fcsgx_metadata,
)

# FCS-GX taxonomy loaders
from .fcsgx_taxonomy_loaders import (
    load_fcsgx_taxonomy_as_table,
)

# FCS-GX report loaders
from .fcsgx_report_loaders import (
    load_fcsgx_results,
    load_fcsgx_report_as_table,
)

# K-mer loaders
from .kmer_loaders import (
    load_kmer_dim_reduction_results,
)

# Reference loaders
from .reference_loaders import (
    process_reference_file_line_by_line,
)

# Utility functions
from .utils import (
    find_files_in_dir,
)

# Make all loaders available at the package level
__all__ = [
    # Basic loaders
    'load_samplesheet',
    'load_yaml_params',
    'load_barcode_check_results',
    
    # Table loaders
    'load_contamination_check_merged_table',
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
