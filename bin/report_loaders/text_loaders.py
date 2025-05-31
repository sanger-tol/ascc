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
    """Load autofiltering results.
    
    This function now uses the table loader to format the autofiltering data as an HTML table.
    """
    # Import the table loader function
    from .table_loaders import load_autofiltering_results_as_table
    
    # Use the table loader function instead
    return load_autofiltering_results_as_table(file_path)
