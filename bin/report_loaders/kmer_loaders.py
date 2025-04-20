#!/usr/bin/env python3
"""
K-mer dimensionality reduction loaders for ASCC HTML report generation.

This module contains functions for loading and formatting k-mer dimensionality reduction results.
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
    if not os.path.exists(dir_path) or not os.path.isdir(dir_path):
        print(f"Directory not found or not a directory: {dir_path}", file=sys.stderr)
        return None
    
    # Check if directory is empty
    if not os.listdir(dir_path):
        print(f"Directory is empty: {dir_path}", file=sys.stderr)
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

            # Look for visualisation files
            viz_files = [f for f in os.listdir(method_dir) if f.endswith("_visualisation.png")]
            for viz_file in viz_files:
                try:
                    with open(os.path.join(method_dir, viz_file), "rb") as f:
                        image_data = f.read()
                        results[f"{method_name}_{viz_file}"] = base64.b64encode(image_data).decode("utf-8")
                except Exception as e:
                    print(f"Error reading visualisation file {viz_file}: {e}", file=sys.stderr)

        return results

    except Exception as e:
        print(f"Error loading kmer dim reduction results: {e}", file=sys.stderr)
        return None
