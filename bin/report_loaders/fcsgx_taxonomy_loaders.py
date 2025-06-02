#!/usr/bin/env python3
"""
FCS-GX taxonomy loaders for ASCC HTML report generation.

This module contains functions for loading and formatting FCS-GX taxonomy reports.
"""

import os
import sys
import pandas as pd
import tempfile


def _fallback_load_fcsgx_taxonomy(file_path, metadata=None):
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
        with open(file_path, "r") as f:
            for i, line in enumerate(f):
                if i >= 100:  # Only read the first 100 lines
                    break
                if not (line.startswith("#") or "GX taxonomy analysis report" in line):
                    sample_lines.append(line.strip())

        if not sample_lines:
            print(
                f"No valid data lines found in taxonomy file: {file_path}",
                file=sys.stderr,
            )
            return metadata, None

        # Create a simple HTML table from the sample
        table_html = "<table class='table table-striped' id='fcsgx_taxonomy_table'>"

        # Add a note that this is a sample
        table_html += "<caption>Note: This is a sample of the first 100 lines due to parsing issues with the full file.</caption>"

        # Add rows
        for line in sample_lines:
            table_html += "<tr>"
            for cell in line.split("\t"):
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


def load_fcsgx_taxonomy_as_table(file_path):
    """Load FCS-GX taxonomy report file and convert to metadata + HTML table.

    Args:
        file_path: Path to the FCS-GX taxonomy report file (*.taxonomy.rpt)

    Returns:
        tuple: (metadata_html, table_html) where metadata_html is the formatted metadata
               and table_html is the HTML table of the tabular data
    """
    if not os.path.exists(file_path):
        print(f"FCS-GX taxonomy file does not exist: {file_path}", file=sys.stderr)
        return None, None
    if os.path.getsize(file_path) == 0:
        print(f"FCS-GX taxonomy file is empty: {file_path}", file=sys.stderr)
        return None, None

    try:
        # Read the first few lines to extract metadata and header
        metadata_raw = None
        header_line = None

        with open(file_path, "r") as f:
            # Read the first 10 lines to find metadata and header
            for i, line in enumerate(f):
                if i >= 10:  # Only check the first 10 lines
                    break

                if "GX taxonomy analysis report" in line:
                    metadata_raw = line.strip()
                    print(
                        f"Found taxonomy metadata at line {i}: {metadata_raw[:50]}...",
                        file=sys.stderr,
                    )
                elif line.startswith("seq-id") or line.startswith("#seq-id"):
                    header_line = line.strip()
                    if header_line.startswith("#"):
                        header_line = header_line[1:]  # Remove leading # if present
                    print(
                        f"Found taxonomy header at line {i}: {header_line}",
                        file=sys.stderr,
                    )
                    break

        # Format the metadata as a table
        metadata_table = None
        if metadata_raw:
            from .fcsgx_loaders import format_fcsgx_metadata

            metadata_table = format_fcsgx_metadata(
                metadata_raw, "fcsgx_taxonomy_metadata_table"
            )
            if metadata_table:
                print(
                    f"Successfully formatted FCS-GX taxonomy metadata as table",
                    file=sys.stderr,
                )
            else:
                print(
                    f"Failed to format FCS-GX taxonomy metadata as table, will use raw text",
                    file=sys.stderr,
                )
                metadata_table = metadata_raw

        if not header_line:
            print(
                f"Could not find header line in taxonomy file: {file_path}",
                file=sys.stderr,
            )
            return metadata_table, None

        # Process the file in chunks to handle large files
        print(f"Processing taxonomy file in chunks: {file_path}", file=sys.stderr)

        # Create a temporary file to store processed data
        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_file:
            # Write the header
            temp_file.write(header_line + "\n")

            # Process the file in chunks, skipping metadata and header lines
            with open(file_path, "r") as f:
                for line in f:
                    # Skip metadata and header lines
                    if (
                        "GX taxonomy analysis report" in line
                        or line.startswith("seq-id")
                        or line.startswith("#seq-id")
                    ):
                        continue

                    # Write data lines to temp file
                    if line.strip():
                        temp_file.write(line)

            temp_file_path = temp_file.name

        # Read the processed data into a DataFrame
        try:
            print(f"Reading processed taxonomy data from temp file", file=sys.stderr)
            df = pd.read_csv(temp_file_path, sep="\t")

            # Remove the temporary file
            os.unlink(temp_file_path)

            print(
                f"Successfully loaded taxonomy data with {len(df)} rows",
                file=sys.stderr,
            )

            # Remove separator columns (sep1, sep2, sep3, sep4, sep5) that only contain pipe symbols
            separator_columns = [col for col in df.columns if col.startswith("sep")]
            if separator_columns:
                print(
                    f"Removing separator columns: {separator_columns}", file=sys.stderr
                )
                df = df.drop(columns=separator_columns)

            # Convert DataFrame to HTML table
            table_html = df.to_html(
                classes="table table-striped",
                index=False,
                table_id="fcsgx_taxonomy_table",
            )

            return metadata_table, table_html

        except Exception as e:
            # Clean up temp file if there was an error
            if os.path.exists(temp_file_path):
                os.unlink(temp_file_path)
            print(f"Error converting taxonomy data to DataFrame: {e}", file=sys.stderr)

            # Fallback to a simpler approach for problematic files
            print(
                f"Falling back to simpler parsing approach for taxonomy file",
                file=sys.stderr,
            )
            return _fallback_load_fcsgx_taxonomy(file_path, metadata_table)

    except Exception as e:
        print(f"Error parsing FCS-GX taxonomy report: {e}", file=sys.stderr)
        return None, None
