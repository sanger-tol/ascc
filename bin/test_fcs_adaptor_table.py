#!/usr/bin/env python3
"""
Test script for the FCS-adaptor table loader.

This script tests the load_fcs_adaptor_results function to ensure it correctly
formats FCS-adaptor data as an HTML table.
"""

import sys
import os

# Add the parent directory to the Python path so we can import the report_loaders package
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import the function to test
from report_loaders import load_fcs_adaptor_results

# Paths to the test FCS-adaptor files
fcs_adaptor_euk_file = "/home/ea/sanger_data/tol_assembly/20250103_ascc/testruns/jinja/run18/test_00/chlamydomonas_dataset_PRIMARY/FCS-adaptor/chlamydomonas_dataset_PRIMARY_euk.fcs_adaptor_report.txt"
fcs_adaptor_prok_file = "/home/ea/sanger_data/tol_assembly/20250103_ascc/testruns/jinja/run18/test_00/chlamydomonas_dataset_PRIMARY/FCS-adaptor/chlamydomonas_dataset_PRIMARY_prok.fcs_adaptor_report.txt"

# Test the function with the eukaryotic file
print("Testing FCS-adaptor table loader with eukaryotic file...")
try:
    # Load and format the FCS-adaptor data
    fcs_adaptor_euk_table = load_fcs_adaptor_results(fcs_adaptor_euk_file)
    
    # Check if the result contains HTML table elements
    if "<table" in fcs_adaptor_euk_table and "<tr" in fcs_adaptor_euk_table and "<td" in fcs_adaptor_euk_table:
        print("✅ Successfully formatted eukaryotic FCS-adaptor data as an HTML table")
        
        # Save the HTML table to a file for inspection
        output_file = "fcs_adaptor_euk_table_test.html"
        with open(output_file, "w") as f:
            f.write("""
            <!DOCTYPE html>
            <html>
            <head>
                <title>FCS-adaptor Eukaryotic Table Test</title>
                <style>
                    body {
                        font-family: sans-serif;
                        max-width: 1200px;
                        margin: 0 auto;
                        padding: 20px;
                    }
                    table {
                        border-collapse: collapse;
                        width: 100%;
                        margin: 1em 0;
                    }
                    th {
                        background-color: #edf2f7;
                        color: #2c5282;
                        font-weight: bold;
                        text-align: left;
                        padding: 8px;
                        border: 1px solid #e2e8f0;
                    }
                    td {
                        padding: 8px;
                        border: 1px solid #e2e8f0;
                    }
                    tr:nth-child(even) {
                        background-color: #f7fafc;
                    }
                    .table-responsive {
                        width: 100%;
                        overflow-x: auto;
                        display: block;
                        margin-bottom: 1em;
                        border: 1px solid #e2e8f0;
                        border-radius: 5px;
                        padding-bottom: 15px;
                    }
                </style>
            </head>
            <body>
                <h1>FCS-adaptor Eukaryotic Table Test</h1>
                """ + fcs_adaptor_euk_table + """
            </body>
            </html>
            """)
        print(f"Saved HTML table to {output_file}")
    else:
        print("❌ Failed to format eukaryotic FCS-adaptor data as an HTML table")
        print(f"Result: {fcs_adaptor_euk_table[:100]}...")
except Exception as e:
    print(f"❌ Error testing eukaryotic FCS-adaptor table loader: {e}")

# Test the function with the prokaryotic file
print("\nTesting FCS-adaptor table loader with prokaryotic file...")
try:
    # Load and format the FCS-adaptor data
    fcs_adaptor_prok_table = load_fcs_adaptor_results(fcs_adaptor_prok_file)
    
    # For the prokaryotic file, we expect a message since it only contains a header
    if "FCS-Adaptor check was run, but no adaptor contamination was detected." in fcs_adaptor_prok_table:
        print("✅ Correctly handled empty prokaryotic FCS-adaptor file")
    else:
        print("❌ Failed to handle empty prokaryotic FCS-adaptor file correctly")
        print(f"Result: {fcs_adaptor_prok_table[:100]}...")
except Exception as e:
    print(f"❌ Error testing prokaryotic FCS-adaptor table loader: {e}")

print("\nTest completed.")
