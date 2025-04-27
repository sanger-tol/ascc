#!/usr/bin/env python3
"""
Fallback functions for parsing FCS-GX taxonomy reports when the main methods fail.
"""

import os
import sys


def fallback_load_fcsgx_taxonomy(file_path, metadata=None):
    """Fallback method to load FCS-GX taxonomy report when the main method fails.
    
    This method uses a simpler approach that may work better for problematic files.
    
    Args:
        file_path: Path to the FCS-GX taxonomy report file
        metadata: The already extracted metadata (if available)
        
    Returns:
        tuple: (metadata_html, table_html) where metadata_html is the formatted metadata
               and table_html is the HTML table of the tabular data
    """
    try:
        # Read the first 100 lines of the file to extract a sample
        sample_lines = []
        with open(file_path, 'r') as f:
            for i, line in enumerate(f):
                if i >= 100:  # Only read the first 100 lines
                    break
                if not (line.startswith("#") or "GX taxonomy analysis report" in line):
                    sample_lines.append(line.strip())
        
        if not sample_lines:
            print(f"No valid data lines found in taxonomy file: {file_path}", file=sys.stderr)
            return metadata, None
            
        # Create a simple HTML table from the sample
        table_html = "<table class='table table-striped' id='fcsgx_taxonomy_table'>"
        
        # Add a note that this is a sample
        table_html += "<caption>Note: This is a sample of the first 100 lines due to parsing issues with the full file.</caption>"
        
        # Add rows
        for line in sample_lines:
            table_html += "<tr>"
            for cell in line.split('\t'):
                table_html += f"<td>{cell}</td>"
            table_html += "</tr>"
        
        table_html += "</table>"
        
        # Wrap table for responsive display
        wrapped_table_html = f"""
        <div class="outer-container">
            <div class="table-responsive">
                <div class="table-wrapper">
                    {table_html}
                </div>
            </div>
        </div>
        """
        
        return metadata, wrapped_table_html
        
    except Exception as e:
        print(f"Error in fallback taxonomy parsing: {e}", file=sys.stderr)
        return metadata, None
