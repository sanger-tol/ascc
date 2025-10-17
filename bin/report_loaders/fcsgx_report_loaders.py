#!/usr/bin/env python3
"""
FCS-GX report loaders for ASCC HTML report generation.

This module contains functions for loading and formatting FCS-GX report data.
Developed by Eerik Aunin (eeaunin@gmail.com)
"""

import os
import sys
import pandas as pd
import tempfile
from .fcsgx_loaders import format_fcsgx_metadata


def load_fcsgx_results(file_path):
    """Load FCS-GX results as plain text.

    Args:
        file_path (str): Path to the FCS-GX report file

    Returns:
        str: Raw content of the FCS-GX report, or an error message if loading fails
    """
    if not os.path.exists(file_path):
        print(f"FCS-GX file does not exist: {file_path}", file=sys.stderr)
        return "FCS-GX check was not run."
    if os.path.getsize(file_path) == 0:
        print(f"FCS-GX file is empty: {file_path}", file=sys.stderr)
        return "FCS-GX check was run and no issues were found."
    try:
        with open(file_path, "r") as f:
            data = f.read()
        # If the file contains only whitespace or just a header line, consider it as "no issues found"
        if not data.strip() or data.strip() == "accessionlengthactionrangename":
            print(
                f"FCS-GX file contains only whitespace or header: {file_path}",
                file=sys.stderr,
            )
            return "FCS-GX check was run and no issues were found."
        # Add debug output
        print(
            f"Successfully loaded FCS-GX data from {file_path}, content length: {len(data)}",
            file=sys.stderr,
        )
        print(f"First 100 characters: {data[:100]}", file=sys.stderr)
        # Return the raw content for display in the report
        return data
    except Exception as e:
        print(f"Error loading FCS-GX results: {e}", file=sys.stderr)
        return "Error loading FCS-GX check results."


def load_fcsgx_report_as_table(file_path):
    """Load FCS-GX report file and convert to metadata + HTML table.

    Args:
        file_path (str): Path to the FCS-GX report file

    Returns:
        tuple: (metadata_html, table_html) where metadata_html is the formatted metadata
               and table_html is the HTML table of the tabular data
    """
    if not os.path.exists(file_path):
        print(f"FCS-GX report file does not exist: {file_path}", file=sys.stderr)
        return None, None
    if os.path.getsize(file_path) == 0:
        print(f"FCS-GX report file is empty: {file_path}", file=sys.stderr)
        return None, None

    try:
        # Read the entire file to extract metadata, header, and data
        metadata_raw = None
        header_line = None
        data_lines = []

        with open(file_path, "r") as f:
            lines = f.readlines()

        # Process each line to identify metadata, header, and data
        for i, line in enumerate(lines):
            line = line.strip()
            if not line:
                continue

            if line.startswith("##"):
                metadata_raw = line
                print(
                    f"Found FCS-GX metadata at line {i}: {metadata_raw[:50]}...",
                    file=sys.stderr,
                )
            elif line.startswith("#"):
                # This is a header line
                header_line = line[1:]  # Remove the leading #
                print(f"Found header at line {i}: {header_line}", file=sys.stderr)
            else:
                # This is a data line
                data_lines.append(line)

        # Format the metadata as a table
        metadata_table = None
        if metadata_raw:
            metadata_table = format_fcsgx_metadata(
                metadata_raw, "fcsgx_report_metadata_table"
            )
            if metadata_table:
                print(
                    f"Successfully formatted FCS-GX report metadata as table",
                    file=sys.stderr,
                )
            else:
                print(
                    f"Failed to format FCS-GX report metadata as table, will use raw text",
                    file=sys.stderr,
                )
                metadata_table = f"<pre><code>{metadata_raw}</code></pre>"

        # Process the tabular data
        table_html = None

        # If we have a header line, try to parse the data as a table
        if header_line:
            try:
                # Create a temporary file with the header and data
                with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_file:
                    temp_file.write(header_line + "\n")
                    for line in data_lines:
                        temp_file.write(line + "\n")

                    temp_file_path = temp_file.name

                # Read the processed data into a DataFrame
                df = pd.read_csv(temp_file_path, sep="\t")

                # Remove the temporary file
                os.unlink(temp_file_path)

                # Check if we have any data rows
                if len(df) > 0:
                    print(
                        f"Successfully loaded FCS-GX report data with {len(df)} rows",
                        file=sys.stderr,
                    )

                    # Remove separator columns (sep1, sep2, sep3, sep4, sep5) that only contain pipe symbols
                    separator_columns = [
                        col for col in df.columns if col.startswith("sep")
                    ]
                    if separator_columns:
                        print(
                            f"Removing separator columns: {separator_columns}",
                            file=sys.stderr,
                        )
                        df = df.drop(columns=separator_columns)

                    # Convert DataFrame to HTML table
                    table_html = df.to_html(
                        classes="table table-striped",
                        index=False,
                        table_id="fcsgx_report_table",
                    )

                    # No longer wrapping the table here, as it will be wrapped in the template
                else:
                    print(f"No data rows found in FCS-GX report", file=sys.stderr)
                    table_html = "<p>No data rows found in the FCS-GX report file.</p>"
            except Exception as e:
                print(
                    f"Error converting FCS-GX report data to DataFrame: {e}",
                    file=sys.stderr,
                )
                # Fall back to raw text for the data section
                if data_lines:
                    data_content = "\n".join(data_lines)
                    table_html = f"<pre><code>{data_content}</code></pre>"
                else:
                    table_html = "<p>No data rows found in the FCS-GX report file.</p>"
        else:
            # If we couldn't identify a header, check if we have any data lines
            if data_lines:
                data_content = "\n".join(data_lines)
                table_html = f"<pre><code>{data_content}</code></pre>"
            else:
                table_html = "<p>No data rows found in the FCS-GX report file.</p>"

        return metadata_table, table_html

    except Exception as e:
        print(f"Error parsing FCS-GX report: {e}", file=sys.stderr)
        return None, None
