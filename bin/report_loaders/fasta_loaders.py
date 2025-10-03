#!/usr/bin/env python3
"""
FASTA loaders for ASCC HTML report generation.

This module contains functions for loading and formatting FASTA-related data,
such as sanitation logs and length filtering logs.
Developed by Eerik Aunin (eeaunin@gmail.com)
"""

import os
import sys
import json
import pandas as pd


def load_fasta_length_filtering_log(file_path):
    """Load FASTA length filtering log data."""
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return None
    try:
        with open(file_path, "r") as f:
            data = json.load(f)

        # Create a summary table from the log data
        summary_data = [
            {
                "Metric": "Total sequences processed",
                "Value": f"{data.get('total_sequences', 0):,}",
            },
            {
                "Metric": "Sequences retained",
                "Value": f"{data.get('sequences_retained', 0):,}",
            },
            {
                "Metric": "Sequences filtered",
                "Value": f"{data.get('sequences_filtered', 0):,}",
            },
            {
                "Metric": "Filter mode",
                "Value": data.get("filter_mode", "unknown").replace("_", " ").title(),
            },
            {"Metric": "Cutoff value", "Value": f"{data.get('cutoff_value', 0):,}"},
        ]

        # Convert to HTML table
        df = pd.DataFrame(summary_data)
        summary_table = df.to_html(
            classes="table table-striped",
            index=False,
            table_id="fasta_length_filtering_summary_table",
        )

        # Create a detailed changes table if there are any changes
        detailed_table = None
        if data.get("detailed_changes"):
            # Convert detailed changes to a DataFrame
            detailed_data = []
            for change in data.get("detailed_changes", []):
                row = {
                    "Header": change.get("header", ""),
                    "Change Type": change.get("change_type", "")
                    .replace("_", " ")
                    .title(),
                }

                if "length" in change:
                    row["Length"] = f"{change['length']:,}"

                if "reason" in change:
                    row["Reason"] = change["reason"]

                if "details" in change:
                    row["Details"] = change["details"]

                detailed_data.append(row)

            # Convert to HTML table
            detailed_df = pd.DataFrame(detailed_data)
            detailed_table = detailed_df.to_html(
                classes="table table-striped",
                index=False,
                table_id="fasta_length_filtering_detailed_table",
            )

        # Create a more compact table without forcing horizontal scrolling
        summary_table_html = f"""
        <div class="fasta-sanitation-table">
            {summary_table}
        </div>
        """

        detailed_table_html = ""
        if detailed_table:
            detailed_table_html = f"""
            <h4>Detailed Filtering</h4>
            <div class="fasta-sanitation-table">
                {detailed_table}
            </div>
            """

        # Combine the tables with a status indicator
        has_filtering = data.get("has_filtering", False)
        status_class = (
            "notice-container notice-error"
            if has_filtering
            else "notice-container notice-success"
        )
        status_message = (
            f"{data.get('sequences_filtered', 0):,} sequences were filtered by length."
            if has_filtering
            else "No sequences were filtered by length."
        )

        result_html = f"""
        <div class="{status_class}">
            <p>{status_message}</p>
        </div>
        <h4>Length Filtering Summary</h4>
        {summary_table_html}
        {detailed_table_html}
        """

        return {"html": result_html, "has_filtering": has_filtering}
    except Exception as e:
        print(f"Error loading FASTA length filtering log: {e}", file=sys.stderr)
        return None


def load_fasta_sanitation_log(file_path, length_filtering_file_path=None):
    """Load FASTA sanitation log data and optionally combine with length filtering data."""
    sanitation_data = None
    length_filtering_data = None

    # Load sanitation log
    if file_path and os.path.exists(file_path) and os.path.getsize(file_path) > 0:
        try:
            with open(file_path, "r") as f:
                data = json.load(f)

            # Create a summary table from the log data
            summary_data = [
                {
                    "Metric": "Total sequences processed",
                    "Value": f"{data.get('total_sequences', 0):,}",
                },
                {
                    "Metric": "Headers shortened",
                    "Value": f"{data.get('headers_shortened', 0):,}",
                },
                {
                    "Metric": "Headers with problematic characters",
                    "Value": f"{data.get('headers_with_problematic_chars', 0):,}",
                },
                {
                    "Metric": "Headers with commas",
                    "Value": f"{data.get('headers_with_commas', 0):,}",
                },
                {
                    "Metric": "Sequences with non-ATGC bases",
                    "Value": f"{data.get('sequences_with_non_atgc', 0):,}",
                },
                {
                    "Metric": "Non-ATGC bases replaced",
                    "Value": f"{data.get('non_atgc_bases_replaced', 0):,}",
                },
                {
                    "Metric": "All-N sequences skipped",
                    "Value": f"{data.get('all_n_sequences_skipped', 0):,}",
                },
                {
                    "Metric": "Duplicate headers detected",
                    "Value": f"{data.get('duplicate_headers_detected', 0):,}",
                },
            ]

            # Convert to HTML table
            df = pd.DataFrame(summary_data)
            summary_table = df.to_html(
                classes="table table-striped",
                index=False,
                table_id="fasta_sanitation_summary_table",
            )

            # Create a detailed changes table if there are any changes
            detailed_table = None
            if data.get("detailed_changes"):
                # Convert detailed changes to a DataFrame
                detailed_data = []
                for change in data.get("detailed_changes", []):
                    row = {
                        "Header": change.get("header", ""),
                        "Change Type": change.get("change_type", "")
                        .replace("_", " ")
                        .title(),
                    }

                    if "original" in change and "modified" in change:
                        row["Original"] = change["original"]
                        row["Modified"] = change["modified"]
                    elif "count" in change:
                        row["Count"] = change["count"]

                    if "details" in change:
                        row["Details"] = change["details"]

                    detailed_data.append(row)

                # Convert to HTML table
                detailed_df = pd.DataFrame(detailed_data)
                detailed_table = detailed_df.to_html(
                    classes="table table-striped",
                    index=False,
                    table_id="fasta_sanitation_detailed_table",
                )

            # Create a more compact table without forcing horisontal scrolling
            summary_table_html = f"""
            <div class="fasta-sanitation-table">
                {summary_table}
            </div>
            """

            detailed_table_html = ""
            if detailed_table:
                detailed_table_html = f"""
                <h4>Detailed Header Changes</h4>
                <div class="fasta-sanitation-table">
                    {detailed_table}
                </div>
                """

            # Combine the tables with a status indicator
            has_issues = data.get("has_issues", False)
            status_class = (
                "notice-container notice-error"
                if has_issues
                else "notice-container notice-success"
            )
            status_message = (
                "Issues with the input FASTA file headers were detected."
                if has_issues
                else "No issues with the input FASTA file headers were detected."
            )

            sanitation_html = f"""
            <div class="{status_class}">
                <p>{status_message}</p>
            </div>
            <h4>Header Sanitation Summary</h4>
            {summary_table_html}
            {detailed_table_html}
            """

            sanitation_data = {"html": sanitation_html, "has_issues": has_issues}
        except Exception as e:
            print(f"Error loading FASTA sanitation log: {e}", file=sys.stderr)

    # Load length filtering log if provided
    if (
        length_filtering_file_path
        and os.path.exists(length_filtering_file_path)
        and os.path.getsize(length_filtering_file_path) > 0
    ):
        length_filtering_data = load_fasta_length_filtering_log(
            length_filtering_file_path
        )

    # Combine the results
    if sanitation_data and length_filtering_data:
        combined_html = f"""
        <h3>FASTA Header Sanitation</h3>
        {sanitation_data["html"]}
        <h3>FASTA Length Filtering</h3>
        {length_filtering_data["html"]}
        """

        return {
            "html": combined_html,
            "has_issues": sanitation_data["has_issues"]
            or length_filtering_data["has_filtering"],
        }
    elif sanitation_data:
        return sanitation_data
    elif length_filtering_data:
        return length_filtering_data
    else:
        return None
