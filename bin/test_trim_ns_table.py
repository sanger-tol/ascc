#!/usr/bin/env python3
"""
Test script for the Trailing Ns table loader.

This script tests the load_trim_Ns_results function to ensure it correctly
formats Trailing Ns data as HTML tables, including the warning comments.
"""

import sys
import os

# Add the parent directory to the Python path so we can import the report_loaders package
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import the function to test
from report_loaders import load_trim_Ns_results

# Path to the test Trailing Ns file
trim_ns_file = "/home/ea/sanger_data/tol_assembly/20250103_ascc/testruns/jinja/run18/work/40/cd77d8c31a32d29376093a3803a081/chlamydomonas_dataset_PRIMARY_trim_Ns"

# Test the function
print("Testing Trailing Ns table loader...")
try:
    # Load and format the Trailing Ns data
    trim_ns_html = load_trim_Ns_results(trim_ns_file)
    
    # Check if the result contains HTML table elements
    if "<table" in trim_ns_html and "<tr" in trim_ns_html and "<td" in trim_ns_html:
        print("✅ Successfully formatted Trailing Ns data as HTML tables")
        
        # Check if both tables are present
        if "trim_ns_warnings_table" in trim_ns_html and "trim_ns_table" in trim_ns_html:
            print("✅ Both warnings table and actions table are present")
        elif "trim_ns_warnings_table" in trim_ns_html:
            print("✅ Warnings table is present")
            print("❌ Actions table is missing")
        elif "trim_ns_table" in trim_ns_html:
            print("❌ Warnings table is missing")
            print("✅ Actions table is present")
        else:
            print("❌ Neither table is present")
        
        # Save the HTML to a file for inspection
        output_file = "trim_ns_table_test.html"
        with open(output_file, "w") as f:
            f.write("""
            <!DOCTYPE html>
            <html>
            <head>
                <title>Trailing Ns Table Test</title>
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
                    h4 {
                        color: #2c5282;
                        margin-top: 1.5em;
                        border-bottom: 1px solid #e2e8f0;
                        padding-bottom: 0.3em;
                    }
                </style>
            </head>
            <body>
                <h1>Trailing Ns Table Test</h1>
                """ + trim_ns_html + """
            </body>
            </html>
            """)
        print(f"Saved HTML to {output_file}")
    else:
        print("❌ Failed to format Trailing Ns data as HTML tables")
        print(f"Result: {trim_ns_html[:100]}...")
except Exception as e:
    print(f"❌ Error testing Trailing Ns table loader: {e}")

print("\nTest completed.")
