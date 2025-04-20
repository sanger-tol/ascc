#!/usr/bin/env python3
"""
Test script for the VecScreen table loader.

This script tests the load_vecscreen_results function to ensure it correctly
formats VecScreen data as an HTML table.
"""

import sys
import os

# Add the parent directory to the Python path so we can import the report_loaders package
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import the function to test
from report_loaders import load_vecscreen_results

# Path to the test VecScreen file
vecscreen_file = "/home/ea/sanger_data/tol_assembly/20250103_ascc/testruns/jinja/run18/test_00/chlamydomonas_dataset_PRIMARY/summarise_vecscreen_output/chlamydomonas_dataset_PRIMARY.vecscreen_contamination"

# Test the function
print("Testing VecScreen table loader...")
try:
    # Load and format the VecScreen data
    vecscreen_table = load_vecscreen_results(vecscreen_file)
    
    # Check if the result contains HTML table elements
    if "<table" in vecscreen_table and "<tr" in vecscreen_table and "<td" in vecscreen_table:
        print("✅ Successfully formatted VecScreen data as an HTML table")
        
        # Save the HTML table to a file for inspection
        output_file = "vecscreen_table_test.html"
        with open(output_file, "w") as f:
            f.write("""
            <!DOCTYPE html>
            <html>
            <head>
                <title>VecScreen Table Test</title>
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
                <h1>VecScreen Table Test</h1>
                """ + vecscreen_table + """
            </body>
            </html>
            """)
        print(f"Saved HTML table to {output_file}")
    else:
        print("❌ Failed to format VecScreen data as an HTML table")
        print(f"Result: {vecscreen_table[:100]}...")
except Exception as e:
    print(f"❌ Error testing VecScreen table loader: {e}")

print("\nTest completed.")
