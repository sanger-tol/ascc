#!/usr/bin/env python3
"""
Reference file loaders for ASCC HTML report generation.

This module contains functions for loading and processing reference genome files.
"""

import os
import sys
import pandas as pd


def process_reference_file_line_by_line(file_path):
    """Process reference file line by line to avoid loading large files into memory.

    Returns a list of dictionaries with assembly statistics that can be converted to an HTML table.

    Args:
        file_path (str): Path to the reference FASTA file

    Returns:
        str: HTML table with assembly statistics, or None if processing fails
    """
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return None

    try:
        # Extract basic information without loading the entire file
        sequence_count = 0
        total_length = 0
        gc_count = 0
        n_count = 0

        with open(file_path, "r") as f:
            current_seq = ""
            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    # Process the previous sequence if it exists
                    if current_seq:
                        seq_length = len(current_seq)
                        total_length += seq_length
                        gc_count += current_seq.upper().count(
                            "G"
                        ) + current_seq.upper().count("C")
                        n_count += current_seq.upper().count("N")

                    sequence_count += 1
                    current_seq = ""
                else:
                    current_seq += line

            # Process the last sequence
            if current_seq:
                seq_length = len(current_seq)
                total_length += seq_length
                gc_count += current_seq.upper().count("G") + current_seq.upper().count(
                    "C"
                )
                n_count += current_seq.upper().count("N")

        # Calculate GC content
        non_n_bases = total_length - n_count
        gc_content = (gc_count / non_n_bases * 100) if non_n_bases > 0 else 0
        n_percent = (n_count / total_length * 100) if total_length > 0 else 0

        # Create a list of dictionaries for the table
        stats_list = [
            {"Statistic": "Sequence count", "Value": f"{sequence_count:,}"},
            {"Statistic": "Total length", "Value": f"{total_length:,} bp"},
            {"Statistic": "GC content", "Value": f"{gc_content:.2f}%"},
            {"Statistic": "N count", "Value": f"{n_count:,} ({n_percent:.2f}%)"},
        ]

        # Convert to HTML table
        df = pd.DataFrame(stats_list)
        table_html = df.to_html(
            classes="table table-striped", index=False, table_id="assembly_stats_table"
        )

        # Wrap the table in multiple container layers for better scrolling
        return f"""
        <div class="outer-container">
            <div class="table-responsive">
                <div class="table-wrapper">
                    {table_html}
                </div>
            </div>
        </div>
        """
    except Exception as e:
        print(f"Error processing reference file: {e}", file=sys.stderr)
        return None
