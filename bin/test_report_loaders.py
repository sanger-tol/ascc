#!/usr/bin/env python3
"""
Test script for the report_loaders package.

This script imports and tests various functions from the report_loaders package
to ensure they are working correctly.
"""

import sys
import os

# Add the parent directory to the Python path so we can import the report_loaders package
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Test importing from the report_loaders package directly
print("Testing imports from report_loaders package...")
try:
    from report_loaders import (
        load_samplesheet,
        load_yaml_params,
        load_barcode_check_results,
        load_contamination_check_merged_table,
        load_trim_Ns_results,
        load_fcs_adaptor_results,
        load_vecscreen_results,
        load_autofiltering_results,
        load_fasta_length_filtering_log,
        load_fasta_sanitation_log,
        parse_fcsgx_metadata,
        flatten_nested_dict,
        metadata_to_html_table,
        format_fcsgx_metadata,
        fallback_load_fcsgx_taxonomy,
        load_fcsgx_taxonomy_as_table,
    )
    print("✅ Successfully imported all functions from report_loaders package")
except ImportError as e:
    print(f"❌ Error importing from report_loaders package: {e}")
    sys.exit(1)

# Test importing from the html_report_loaders module (backward compatibility)
print("\nTesting imports from html_report_loaders module...")
try:
    from html_report_loaders import (
        load_samplesheet,
        load_yaml_params,
        load_barcode_check_results,
        load_contamination_check_merged_table,
        load_trim_Ns_results,
        load_fcs_adaptor_results,
        load_vecscreen_results,
        load_autofiltering_results,
        load_fasta_length_filtering_log,
        load_fasta_sanitation_log,
        parse_fcsgx_metadata,
        flatten_nested_dict,
        metadata_to_html_table,
        format_fcsgx_metadata,
        fallback_load_fcsgx_taxonomy,
        load_fcsgx_taxonomy_as_table,
    )
    print("✅ Successfully imported all functions from html_report_loaders module")
except ImportError as e:
    print(f"❌ Error importing from html_report_loaders module: {e}")
    sys.exit(1)

# Test the FCS-GX metadata parsing function
print("\nTesting FCS-GX metadata parsing...")
test_metadata = '##[["FCS genome report", 2, 1], {"git-rev": "v0.5.4-5-g5dfd516", "run-date": "Sat Apr 12 06:02:08 2025"}]'
try:
    parsed = parse_fcsgx_metadata(test_metadata)
    if parsed and "report_info" in parsed and "metadata" in parsed:
        print(f"✅ Successfully parsed FCS-GX metadata: {parsed}")
    else:
        print(f"❌ Failed to parse FCS-GX metadata correctly: {parsed}")
except Exception as e:
    print(f"❌ Error parsing FCS-GX metadata: {e}")

print("\nAll tests completed.")
