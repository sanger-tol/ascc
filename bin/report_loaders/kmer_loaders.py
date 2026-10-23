#!/usr/bin/env python3
"""
K-mer dimensionality reduction loaders for ASCC HTML report generation.

This module contains functions for loading and formatting k-mer dimensionality reduction results.
Developed by Eerik Aunin (eeaunin@gmail.com)
"""

import os
import sys
import base64
import pandas as pd


def load_kmer_dim_reduction_results(dir_path):
    """Load kmer dimensionality reduction results.

    Args:
        dir_path (str): Path to the directory containing k-mer dimensionality reduction results

    Returns:
        dict: Dictionary of results, or None if loading fails
    """
    print(
        f"Loading k-mer dimensionality reduction results from: {dir_path}",
        file=sys.stderr,
    )

    if not os.path.exists(dir_path) or not os.path.isdir(dir_path):
        print(f"Directory not found or not a directory: {dir_path}", file=sys.stderr)
        return None

    # Check if directory is empty
    dir_contents = os.listdir(dir_path)
    if not dir_contents:
        print(f"Directory is empty: {dir_path}", file=sys.stderr)
        return None

    print(f"Directory contents: {dir_contents}", file=sys.stderr)

    results = {}

    # List all directories in the kmers directory
    try:
        # Check if we're looking at a directory containing method directories (or symlinks to them)
        method_dirs = []
        for item in os.listdir(dir_path):
            item_path = os.path.join(dir_path, item)
            target_path = item_path
            is_link_to_dir = False

            # Check if it's a symlink and resolve it
            if os.path.islink(item_path):
                try:
                    resolved_path = os.path.realpath(item_path)  # Follow all symlinks
                    if os.path.isdir(resolved_path):
                        target_path = resolved_path
                        is_link_to_dir = True
                except OSError as e:
                    print(
                        f"Warning: Could not resolve symlink {item_path}: {e}",
                        file=sys.stderr,
                    )
                    continue  # Skip broken links

            # Check if the item itself is a directory OR if it's a link to a directory
            if os.path.isdir(item_path) or is_link_to_dir:
                # Check the name of the target directory (resolved link or the item itself)
                target_basename = os.path.basename(target_path)
                if "_kmers_dim_reduction_dir" in target_basename:
                    # Add the path within the work dir (item_path) to the list
                    method_dirs.append(item_path)

        # If no method directories/links found, assume we're already in a method directory
        if not method_dirs:
            # Check if the current dir_path itself matches the pattern
            dir_basename = os.path.basename(
                os.path.realpath(dir_path)
            )  # Resolve potential link
            if "_kmers_dim_reduction_dir" in dir_basename:
                method_dirs = [dir_path]
            else:
                print(
                    f"No directories matching '*_kmers_dim_reduction_dir' found in {dir_path}",
                    file=sys.stderr,
                )
                return None  # Explicitly return None if no valid dirs found

        print(
            f"Found {len(method_dirs)} dimensionality reduction method sources (dirs/links)",
            file=sys.stderr,
        )

        for (
            method_source_path
        ) in method_dirs:  # Iterate through the paths found (could be dirs or links)
            # Resolve the path to get the actual directory for reading files
            try:
                method_dir_actual = os.path.realpath(method_source_path)
                if not os.path.isdir(method_dir_actual):
                    print(
                        f"Error: Resolved path {method_dir_actual} is not a directory. Skipping {method_source_path}.",
                        file=sys.stderr,
                    )
                    continue
            except OSError as e:
                print(
                    f"Error: Could not resolve path {method_source_path}: {e}. Skipping.",
                    file=sys.stderr,
                )
                continue

            # Extract method name from the basename of the *actual* directory
            dir_basename = os.path.basename(method_dir_actual)
            if "_kmers_dim_reduction_dir" in dir_basename:
                method_parts = dir_basename.replace(
                    "_kmers_dim_reduction_dir", ""
                ).split("_")
                method_name = method_parts[-1]
            else:
                # Fallback if the pattern isn't found in the resolved name (shouldn't happen with the check above)
                method_name = dir_basename.split("_")[0]

            if not method_name or method_name == "kmers":
                method_name = "unknown"  # Assign a default if extraction fails

            print(
                f"Processing method: {method_name} from source {method_source_path} (actual: {method_dir_actual})",
                file=sys.stderr,
            )

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
                            metrics_list.append(
                                {"Metric": metric_name, "Value": metric_value}
                            )

                if not metrics_list:
                    return metrics_content  # Return original content if parsing fails

                # Convert to a pandas DataFrame for display
                df = pd.DataFrame(metrics_list)
                table_html = df.to_html(
                    classes="table table-striped", index=False, table_id=table_id
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

            # Look for metrics files in the *actual* directory
            metrics_files = [
                f
                for f in os.listdir(method_dir_actual)
                if f.endswith("_metrics.txt")
                and not f.endswith("_clustering_metrics.txt")
            ]
            for metrics_file in metrics_files:
                try:
                    with open(os.path.join(method_dir_actual, metrics_file), "r") as f:
                        metrics_content = f.read()
                        # Track which metrics we've loaded
                        for line in metrics_content.splitlines():
                            if ":" in line:
                                metric_name = line.split(":")[0].strip()
                                loaded_metrics.add(metric_name)

                        # Convert to HTML table
                        table_id = f"{method_name}_{os.path.basename(metrics_file).replace('.txt', '')}_table"
                        results[f"{method_name}_{metrics_file}"] = (
                            metrics_to_html_table(metrics_content, table_id)
                        )
                except Exception as e:
                    print(
                        f"Error reading metrics file {metrics_file}: {e}",
                        file=sys.stderr,
                    )

            # Look for clustering metrics files in the *actual* directory
            clustering_files = [
                f
                for f in os.listdir(method_dir_actual)
                if f.endswith("_clustering_metrics.txt")
            ]
            for clustering_file in clustering_files:
                try:
                    with open(
                        os.path.join(method_dir_actual, clustering_file), "r"
                    ) as f:
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
                            results[f"{method_name}_{clustering_file}"] = (
                                metrics_to_html_table("\n".join(new_content), table_id)
                            )
                except Exception as e:
                    print(
                        f"Error reading clustering file {clustering_file}: {e}",
                        file=sys.stderr,
                    )

            # Look for visualisation files in the *actual* directory
            viz_files = [
                f
                for f in os.listdir(method_dir_actual)
                if f.endswith("_visualisation.png")
            ]
            for viz_file in viz_files:
                try:
                    with open(os.path.join(method_dir_actual, viz_file), "rb") as f:
                        image_data = f.read()
                        results[f"{method_name}_{viz_file}"] = base64.b64encode(
                            image_data
                        ).decode("utf-8")
                except Exception as e:
                    print(
                        f"Error reading visualisation file {viz_file}: {e}",
                        file=sys.stderr,
                    )

        return results

    except Exception as e:
        print(f"Error loading kmer dim reduction results: {e}", file=sys.stderr)
        return None
