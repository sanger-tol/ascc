#!/usr/bin/env python3
"""
Module for parsing and formatting FCS-GX report data.
This module contains functions to parse the JSON-like metadata in FCS-GX reports
and convert it to formatted HTML tables.
"""

import json
import re
import pandas as pd
import sys


def parse_fcsgx_metadata(metadata_str):
    """Parse the FCS-GX metadata string into a structured dictionary.
    
    Args:
        metadata_str: Raw metadata string from FCS-GX report
        
    Returns:
        dict: Parsed metadata as a dictionary
    """
    if not metadata_str or not metadata_str.strip():
        return None
        
    # Remove leading '##' if present
    if metadata_str.startswith('##'):
        metadata_str = metadata_str[2:]
        
    try:
        # Try to parse as JSON
        data = json.loads(metadata_str)
        
        # If successful, extract the components
        if isinstance(data, list) and len(data) >= 2:
            # Extract report info and metadata
            report_info = data[0]  # Usually ["FCS genome report", version, subversion]
            metadata = data[1]     # The JSON object with all the metadata
            
            # Create a result dictionary with both parts
            result = {
                "report_info": report_info,
                "metadata": metadata
            }
            return result
    except json.JSONDecodeError:
        # If direct JSON parsing fails, try with regex
        try:
            # Match the pattern: [["report name", version, subversion], {json_data}]
            pattern = r'\[\s*(\[.*?\])\s*,\s*(\{.*\})\s*\]'
            match = re.search(pattern, metadata_str)
            
            if match:
                report_info_str = match.group(1)
                metadata_str = match.group(2)
                
                # Parse each part
                report_info = json.loads(report_info_str)
                metadata = json.loads(metadata_str)
                
                result = {
                    "report_info": report_info,
                    "metadata": metadata
                }
                return result
        except Exception as e:
            print(f"Error parsing FCS-GX metadata with regex: {e}", file=sys.stderr)
    
    # If all parsing attempts fail, return None
    print(f"Failed to parse FCS-GX metadata: {metadata_str[:100]}...", file=sys.stderr)
    return None


def flatten_nested_dict(d, parent_key='', sep='.'):
    """Flatten a nested dictionary into a single level dictionary with dot notation keys.
    
    Args:
        d: The dictionary to flatten
        parent_key: The parent key prefix (used in recursion)
        sep: Separator between nested keys
        
    Returns:
        dict: Flattened dictionary
    """
    items = []
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        
        if isinstance(v, dict):
            items.extend(flatten_nested_dict(v, new_key, sep=sep).items())
        elif isinstance(v, list):
            # Convert lists to comma-separated strings
            items.append((new_key, ", ".join(str(x) for x in v)))
        else:
            items.append((new_key, v))
            
    return dict(items)


def metadata_to_html_table(parsed_metadata, table_id="fcsgx_metadata_table"):
    """Convert parsed metadata to an HTML table.
    
    Args:
        parsed_metadata: The parsed metadata dictionary from parse_fcsgx_metadata()
        table_id: ID to use for the HTML table
        
    Returns:
        str: HTML table representation of the metadata
    """
    if not parsed_metadata:
        return None
        
    table_rows = []
    
    # Add report info rows
    if "report_info" in parsed_metadata and isinstance(parsed_metadata["report_info"], list):
        report_info = parsed_metadata["report_info"]
        if len(report_info) >= 1:
            table_rows.append({"Property": "Report Type", "Value": report_info[0]})
        if len(report_info) >= 2:
            table_rows.append({"Property": "Version", "Value": f"{report_info[1]}"})
        if len(report_info) >= 3:
            table_rows.append({"Property": "Subversion", "Value": f"{report_info[2]}"})
    
    # Add metadata rows
    if "metadata" in parsed_metadata and isinstance(parsed_metadata["metadata"], dict):
        # Flatten the nested metadata
        flat_metadata = flatten_nested_dict(parsed_metadata["metadata"])
        
        # Add each key-value pair to the table
        for key, value in flat_metadata.items():
            table_rows.append({"Property": key, "Value": value})
    
    # Convert to DataFrame and then to HTML
    if table_rows:
        df = pd.DataFrame(table_rows)
        table_html = df.to_html(classes="table table-striped", index=False, table_id=table_id)
        
        # Return the table HTML directly without extra wrappers
        return table_html
    
    return None


def format_fcsgx_metadata(metadata_str, table_id="fcsgx_metadata_table"):
    """Parse and format FCS-GX metadata as an HTML table.
    
    This is the main function to call from other modules.
    
    Args:
        metadata_str: Raw metadata string from FCS-GX report
        table_id: ID to use for the HTML table
        
    Returns:
        str: HTML table representation of the metadata, or None if parsing fails
    """
    parsed_metadata = parse_fcsgx_metadata(metadata_str)
    if parsed_metadata:
        return metadata_to_html_table(parsed_metadata, table_id)
    return None
