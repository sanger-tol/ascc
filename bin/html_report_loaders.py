#!/usr/bin/env python3

import os
import sys
import pandas as pd
import json
import base64
import re


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
        html_table = f"""
        <div class="outer-container">
            <div class="table-responsive">
                <div class="table-wrapper">
                    {table_html}
                </div>
            </div>
        </div>
        """

        return html_table, sample_name
    except Exception as e:
        print(f"Error loading samplesheet: {e}", file=sys.stderr)
        return None, None


def load_yaml_params(file_path=None, params_json=None):
    """Load and format the input parameters as a table.

    This function can load parameters from either:
    1. A YAML file path
    2. A JSON string containing the params object
    """
    params_list = []

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
            return None

        if os.path.getsize(file_path) == 0:
            return None

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
        except Exception as e:
            print(f"Error loading YAML parameters: {e}", file=sys.stderr)

    # If no parameters were found from either source, return a message
    if not params_list:
        return "<p>No parameters found.</p>"

    # Convert to a pandas DataFrame for display
    df = pd.DataFrame(params_list)
    table_html = df.to_html(classes="table table-striped", index=False, table_id="yaml_params_table")

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


def load_barcode_check_results(file_path):
    """Load PacBio barcode check results."""
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return None
    try:
        with open(file_path, "r") as f:
            data = f.read()
        return data
    except Exception as e:
        print(f"Error loading barcode check results: {e}", file=sys.stderr)
        return None


def load_contamination_check_merged_table(file_path):
    """Load contamination check merged table."""
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return None
    try:
        df = pd.read_csv(file_path, sep=",")
        # Add column width styles to make the table wider
        table_html = df.to_html(classes="table table-striped", index=False, table_id="cobiont_check_merged_table")

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
        print(f"Error loading contamination check merged table: {e}", file=sys.stderr)
        return None


def load_fcs_adaptor_results(file_path):
    """Load FCS adaptor results."""
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return None
    try:
        with open(file_path, "r") as f:
            data = f.read()
        return data
    except Exception as e:
        print(f"Error loading FCS adaptor results: {e}", file=sys.stderr)
        return None


def load_trim_Ns_results(file_path):
    """Load trim Ns results."""
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return None
    try:
        with open(file_path, "r") as f:
            data = f.read()
        return data
    except Exception as e:
        print(f"Error loading trim Ns results: {e}", file=sys.stderr)
        return None


def load_vecscreen_results(file_path):
    """Load vecscreen results."""
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return None
    try:
        with open(file_path, "r") as f:
            data = f.read()
        return data
    except Exception as e:
        print(f"Error loading vecscreen results: {e}", file=sys.stderr)
        return None


def load_autofiltering_results(file_path):
    """Load autofiltering results."""
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return None
    try:
        with open(file_path, "r") as f:
            data = f.read()
        return data
    except Exception as e:
        print(f"Error loading autofiltering results: {e}", file=sys.stderr)
        return None


def load_fasta_length_filtering_log(file_path):
    """Load FASTA length filtering log data."""
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return None
    try:
        with open(file_path, "r") as f:
            data = json.load(f)

        # Create a summary table from the log data
        summary_data = [
            {"Metric": "Total sequences processed", "Value": f"{data.get('total_sequences', 0):,}"},
            {"Metric": "Sequences retained", "Value": f"{data.get('sequences_retained', 0):,}"},
            {"Metric": "Sequences filtered", "Value": f"{data.get('sequences_filtered', 0):,}"},
            {"Metric": "Filter mode", "Value": data.get("filter_mode", "unknown").replace("_", " ").title()},
            {"Metric": "Cutoff value", "Value": f"{data.get('cutoff_value', 0):,}"},
        ]

        # Convert to HTML table
        df = pd.DataFrame(summary_data)
        summary_table = df.to_html(
            classes="table table-striped", index=False, table_id="fasta_length_filtering_summary_table"
        )

        # Create a detailed changes table if there are any changes
        detailed_table = None
        if data.get("detailed_changes"):
            # Convert detailed changes to a DataFrame
            detailed_data = []
            for change in data.get("detailed_changes", []):
                row = {
                    "Header": change.get("header", ""),
                    "Change Type": change.get("change_type", "").replace("_", " ").title(),
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
                classes="table table-striped", index=False, table_id="fasta_length_filtering_detailed_table"
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
        status_class = "alert-info" if has_filtering else "alert-success"
        status_message = (
            f"{data.get('sequences_filtered', 0):,} sequences were filtered by length"
            if has_filtering
            else "No sequences were filtered by length"
        )

        result_html = f"""
        <div class="{status_class}">
            {status_message}
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
                {"Metric": "Total sequences processed", "Value": f"{data.get('total_sequences', 0):,}"},
                {"Metric": "Headers shortened", "Value": f"{data.get('headers_shortened', 0):,}"},
                {
                    "Metric": "Headers with problematic characters",
                    "Value": f"{data.get('headers_with_problematic_chars', 0):,}",
                },
                {"Metric": "Headers with commas", "Value": f"{data.get('headers_with_commas', 0):,}"},
                {"Metric": "Sequences with non-ATGC bases", "Value": f"{data.get('sequences_with_non_atgc', 0):,}"},
                {"Metric": "Non-ATGC bases replaced", "Value": f"{data.get('non_atgc_bases_replaced', 0):,}"},
                {"Metric": "All-N sequences skipped", "Value": f"{data.get('all_n_sequences_skipped', 0):,}"},
                {"Metric": "Duplicate headers detected", "Value": f"{data.get('duplicate_headers_detected', 0):,}"},
            ]

            # Convert to HTML table
            df = pd.DataFrame(summary_data)
            summary_table = df.to_html(
                classes="table table-striped", index=False, table_id="fasta_sanitation_summary_table"
            )

            # Create a detailed changes table if there are any changes
            detailed_table = None
            if data.get("detailed_changes"):
                # Convert detailed changes to a DataFrame
                detailed_data = []
                for change in data.get("detailed_changes", []):
                    row = {
                        "Header": change.get("header", ""),
                        "Change Type": change.get("change_type", "").replace("_", " ").title(),
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
                    classes="table table-striped", index=False, table_id="fasta_sanitation_detailed_table"
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
                <h4>Detailed Header Changes</h4>
                <div class="fasta-sanitation-table">
                    {detailed_table}
                </div>
                """

            # Combine the tables with a status indicator
            has_issues = data.get("has_issues", False)
            status_class = "alert-danger" if has_issues else "alert-success"
            status_message = (
                "Issues with the input FASTA file headers were detected"
                if has_issues
                else "No issues with the input FASTA file headers were detected"
            )

            sanitation_html = f"""
            <div class="{status_class}">
                {status_message}
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
        length_filtering_data = load_fasta_length_filtering_log(length_filtering_file_path)

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
            "has_issues": sanitation_data["has_issues"] or length_filtering_data["has_filtering"],
        }
    elif sanitation_data:
        return sanitation_data
    elif length_filtering_data:
        return length_filtering_data
    else:
        return None


def load_fcsgx_results(file_path):
    """Load FCS-GX results."""
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return None
    try:
        with open(file_path, "r") as f:
            data = f.read()
        return data
    except Exception as e:
        print(f"Error loading FCS-GX results: {e}", file=sys.stderr)
        return None


def load_kmer_dim_reduction_results(dir_path):
    """Load kmer dimensionality reduction results."""
    if not os.path.exists(dir_path) or not os.path.isdir(dir_path):
        print(f"Directory not found: {dir_path}", file=sys.stderr)
        return None

    results = {}

    # List all directories in the kmers directory
    try:
        # Check if we're looking at a directory containing method directories
        method_dirs = []
        for item in os.listdir(dir_path):
            item_path = os.path.join(dir_path, item)
            if os.path.isdir(item_path) and "_kmers_dim_reduction_dir" in item:
                method_dirs.append(item_path)

        # If no method directories found, assume we're already in a method directory
        if not method_dirs:
            method_dirs = [dir_path]

        print(f"Found {len(method_dirs)} dimensionality reduction method directories", file=sys.stderr)

        for method_dir in method_dirs:
            # Extract method name by removing the known suffix and taking the last part
            dir_basename = os.path.basename(method_dir)
            if "_kmers_dim_reduction_dir" in dir_basename:
                # Remove the suffix and split by underscores
                method_parts = dir_basename.replace("_kmers_dim_reduction_dir", "").split("_")
                # Take the last element as the method name
                method_name = method_parts[-1]
            else:
                # Fallback to the old method if the expected pattern isn't found
                method_name = dir_basename.split("_")[0]

            if not method_name or method_name == "kmers":
                method_name = "unknown"

            print(f"Processing method: {method_name} in {method_dir}", file=sys.stderr)

            # Track which metrics we've already loaded to avoid duplication
            loaded_metrics = set()

            # First, try to load the combined metrics (if they exist)
            # Function to convert metrics text to HTML table
            def metrics_to_html_table(metrics_content, table_id):
                metrics_list = []
                for line in metrics_content.splitlines():
                    if ":" in line:
                        parts = line.split(":", 1)
                        if len(parts) == 2:
                            metric_name = parts[0].strip()
                            metric_value = parts[1].strip()
                            metrics_list.append({"Metric": metric_name, "Value": metric_value})

                if not metrics_list:
                    return metrics_content  # Return original content if parsing fails

                # Convert to a pandas DataFrame for display
                df = pd.DataFrame(metrics_list)
                table_html = df.to_html(classes="table table-striped", index=False, table_id=table_id)

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

            # Look for metrics files (trustworthiness and continuity)
            metrics_files = [
                f
                for f in os.listdir(method_dir)
                if f.endswith("_metrics.txt") and not f.endswith("_clustering_metrics.txt")
            ]
            for metrics_file in metrics_files:
                try:
                    with open(os.path.join(method_dir, metrics_file), "r") as f:
                        metrics_content = f.read()
                        # Track which metrics we've loaded
                        for line in metrics_content.splitlines():
                            if ":" in line:
                                metric_name = line.split(":")[0].strip()
                                loaded_metrics.add(metric_name)

                        # Convert to HTML table
                        table_id = f"{method_name}_{os.path.basename(metrics_file).replace('.txt', '')}_table"
                        results[f"{method_name}_{metrics_file}"] = metrics_to_html_table(metrics_content, table_id)
                except Exception as e:
                    print(f"Error reading metrics file {metrics_file}: {e}", file=sys.stderr)

            # Look for clustering metrics files (silhouette, calinski_harabasz, etc.)
            clustering_files = [f for f in os.listdir(method_dir) if f.endswith("_clustering_metrics.txt")]
            for clustering_file in clustering_files:
                try:
                    with open(os.path.join(method_dir, clustering_file), "r") as f:
                        clustering_content = f.read()
                        # Only include metrics that haven't been loaded yet
                        new_content = []
                        for line in clustering_content.splitlines():
                            if ":" in line:
                                metric_name = line.split(":")[0].strip()
                                if metric_name not in loaded_metrics:
                                    new_content.append(line)
                                    loaded_metrics.add(metric_name)

                        if new_content:
                            # Convert to HTML table
                            table_id = f"{method_name}_{os.path.basename(clustering_file).replace('.txt', '')}_table"
                            results[f"{method_name}_{clustering_file}"] = metrics_to_html_table(
                                "\n".join(new_content), table_id
                            )
                except Exception as e:
                    print(f"Error reading clustering file {clustering_file}: {e}", file=sys.stderr)

            # Look for visualization files
            viz_files = [f for f in os.listdir(method_dir) if f.endswith("_visualisation.png")]
            for viz_file in viz_files:
                try:
                    with open(os.path.join(method_dir, viz_file), "rb") as f:
                        image_data = f.read()
                        results[f"{method_name}_{viz_file}"] = base64.b64encode(image_data).decode("utf-8")
                except Exception as e:
                    print(f"Error reading visualization file {viz_file}: {e}", file=sys.stderr)

        return results

    except Exception as e:
        print(f"Error loading kmer dim reduction results: {e}", file=sys.stderr)
        return None


def process_reference_file_line_by_line(file_path):
    """Process reference file line by line to avoid loading large files into memory.

    Returns a list of dictionaries with assembly statistics that can be converted to an HTML table.
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
                        gc_count += current_seq.upper().count("G") + current_seq.upper().count("C")
                        n_count += current_seq.upper().count("N")

                    sequence_count += 1
                    current_seq = ""
                else:
                    current_seq += line

            # Process the last sequence
            if current_seq:
                seq_length = len(current_seq)
                total_length += seq_length
                gc_count += current_seq.upper().count("G") + current_seq.upper().count("C")
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
        table_html = df.to_html(classes="table table-striped", index=False, table_id="assembly_stats_table")

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


# Helper function to find files in directories
def find_files_in_dir(directory, pattern=None, extension=None):
    """Find files in a directory matching a pattern or extension."""
    if not directory or not os.path.exists(directory) or not os.path.isdir(directory):
        return []

    files = []
    for f in os.listdir(directory):
        file_path = os.path.join(directory, f)
        if os.path.isfile(file_path):
            if extension and f.endswith(extension):
                files.append(file_path)
            elif pattern and pattern in f:
                files.append(file_path)
            elif not extension and not pattern:
                files.append(file_path)

    return files
