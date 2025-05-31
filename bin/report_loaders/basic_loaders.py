#!/usr/bin/env python3
"""
Basic loaders for ASCC HTML report generation.

This module contains functions for loading and formatting simple data types
like samplesheets, YAML parameters, and other basic file formats.
"""

import os
import sys
import pandas as pd
import json
from .utils import wrap_table_html

# Constants for common messages
NO_BARCODE_CONTAMINATION_MSG = "The check for retained PacBio multiplexing barcode contamination was run. No contamination was found."


def load_samplesheet(file_path):
    """Load and format the input samplesheet CSV.

    Returns a tuple of (html_table, sample_name) where sample_name is extracted from the samplesheet
    if available, otherwise None.
    """
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return None, None
    try:
        df = pd.read_csv(file_path)

        # Try to extract sample name from the samplesheet
        sample_name = None
        if "sample" in df.columns:
            # Get the first sample name
            sample_values = df["sample"].dropna().unique()
            if len(sample_values) > 0:
                sample_name = sample_values[0]

        # Add column width styles to make the table wider
        table_html = df.to_html(classes="table table-striped", index=False, table_id="samplesheet_table")

        # Wrap the table in multiple container layers for better scrolling
        html_table = wrap_table_html(table_html)

        return html_table, sample_name
    except Exception as e:
        print(f"Error loading samplesheet: {e}", file=sys.stderr)
        return None, None


def load_yaml_params(file_path=None, params_json=None):
    """Load and format the input parameters as a table.

    This function can load parameters from either:
    1. A YAML file path
    2. A JSON string containing the params object
    
    Returns:
        tuple: (html_table, params_dict) where html_table is the HTML representation of the parameters
               and params_dict is a dictionary of parameter names to values
    """
    params_list = []
    params_dict = {}  # Dictionary to store parameter values

    # Try to load from JSON string first (this is the new preferred method)
    if params_json:
        try:
            params_dict = json.loads(params_json)

            # Convert the params dictionary to a list of key-value pairs
            for key, value in params_dict.items():
                # Handle nested structures by converting to string
                if isinstance(value, (dict, list)):
                    value_str = json.dumps(value)
                else:
                    value_str = str(value)

                params_list.append({"Parameter": key, "Value": value_str})
        except Exception as e:
            print(f"Error parsing params JSON: {e}", file=sys.stderr)

    # If no params were loaded from JSON, try the YAML file
    if not params_list and file_path:
        if not os.path.exists(file_path):
            return None, {}

        if os.path.getsize(file_path) == 0:
            return None, {}

        try:
            # Read the YAML file as plain text
            with open(file_path, "r") as f:
                yaml_text = f.read()

            # Parse the text to extract key-value pairs
            for line in yaml_text.splitlines():
                line = line.strip()
                # Skip empty lines and comments
                if not line or line.startswith("#"):
                    continue

                # Handle lines with colons (key-value pairs)
                if ":" in line:
                    # Split at the first colon to separate key and value
                    parts = line.split(":", 1)
                    if len(parts) == 2:
                        key = parts[0].strip()
                        value = parts[1].strip()
                        params_list.append({"Parameter": key, "Value": value})
                        params_dict[key] = value  # Store in dictionary
        except Exception as e:
            print(f"Error loading YAML parameters: {e}", file=sys.stderr)

    # If no parameters were found from either source, return a message
    if not params_list:
        return "<p>No parameters found.</p>", {}

    # Convert to a pandas DataFrame for display
    df = pd.DataFrame(params_list)
    table_html = df.to_html(classes="table table-striped", index=False, table_id="yaml_params_table")

    # Wrap the table in multiple container layers for better scrolling
    html_table = wrap_table_html(table_html)
    
    return html_table, params_dict


def load_barcode_check_results(file_path):
    """Load PacBio barcode check results."""
    if not os.path.exists(file_path):
        return "PacBio barcode check was not run."
    if os.path.getsize(file_path) == 0:
        return NO_BARCODE_CONTAMINATION_MSG
    try:
        with open(file_path, "r") as f:
            data = f.read()
        # If the file is empty or contains only whitespace, return the "No data" message
        if not data.strip():
            return NO_BARCODE_CONTAMINATION_MSG
        return data
    except Exception as e:
        print(f"Error loading barcode check results: {e}", file=sys.stderr)
        return "Error loading PacBio barcode check results."
