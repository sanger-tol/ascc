#!/usr/bin/env python3
"""
Text loaders for ASCC HTML report generation.

This module contains functions for loading and formatting text-based data
that doesn't fit into a tabular structure.
"""

import os
import sys


def load_fcs_adaptor_results(file_path):
    """Load FCS adaptor results.
    
    This function now uses the table loader to format the FCS-Adaptor data as an HTML table.
    """
    # Import the table loader function
    from .table_loaders import load_fcs_adaptor_results_as_table
    
    # Use the table loader function instead
    return load_fcs_adaptor_results_as_table(file_path)


def load_vecscreen_results(file_path):
    """Load vecscreen results.
    
    This function now uses the table loader to format the VecScreen data as an HTML table.
    """
    # Import the table loader function
    from .table_loaders import load_vecscreen_results_as_table
    
    # Use the table loader function instead
    return load_vecscreen_results_as_table(file_path)


def load_autofiltering_results(file_path):
    """Load autofiltering results."""
    if not os.path.exists(file_path):
        print(f"Autofiltering file does not exist: {file_path}", file=sys.stderr)
        return "Autofiltering check was not run."
    if os.path.getsize(file_path) == 0:
        print(f"Autofiltering file is empty: {file_path}", file=sys.stderr)
        return "Autofiltering check was run, but no filtering was performed."
    try:
        with open(file_path, "r") as f:
            data = f.read()
        # If the file is empty or contains only whitespace, return the "No data" message
        if not data.strip():
            print(f"Autofiltering file contains only whitespace: {file_path}", file=sys.stderr)
            return "Autofiltering check was run, but no filtering was performed."
        # Add debug output
        print(f"Successfully loaded Autofiltering data from {file_path}, content length: {len(data)}", file=sys.stderr)
        print(f"First 100 characters: {data[:100]}", file=sys.stderr)
        return data
    except Exception as e:
        print(f"Error loading autofiltering results: {e}", file=sys.stderr)
        return "Error loading autofiltering check results."
