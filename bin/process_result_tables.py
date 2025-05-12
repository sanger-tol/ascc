#!/usr/bin/env python3

"""
Script processing the cobiont check result tables to add a combined classification ('merged_classif') column that is based
    on the output of multiple tools. Also generates a table for estimated coverages per 'merged_classif' column
"""

import pandas as pd
import numpy as np
import argparse
import sys
import os


def generate_counts_df(df, counts_df_output_path):
    """
    Creates a table that shows the number of sequences, mean coverage and sequence span per phylum,
    including the source of the classification
    """
    # Group by both merged_classif and merged_classif_source
    grouped = df.groupby(['merged_classif', 'merged_classif_source'])
    
    # Initialize lists to store results
    phylum_list = []
    source_list = []
    count_list = []
    cov_list = []
    span_list = []
    
    # Calculate statistics for each group
    for (classif, source), group in grouped:
        phylum_list.append(classif)
        source_list.append(source)
        count_list.append(len(group))
        cov_values = group["coverage"]
        length_values = group["length"]
        cov_list.append(round(np.mean(cov_values), 2))
        span_sum = sum(length_values) / 1000000
        span_list.append(round(span_sum, 2))
    
    # Create DataFrame with results
    counts_df = pd.DataFrame({
        "phylum": phylum_list,
        "classification_source": source_list,
        "number_of_sequences": count_list,
        "mean_coverage": cov_list,
        "span_mb": span_list
    })
    
    # Save to CSV
    counts_df.to_csv(counts_df_output_path, index=False)


def main(output_folder, sample_id):
    main_results_table_path = "{}/{}_contamination_check_merged_table.csv".format(
        output_folder, sample_id
    )
    btk_results_table_path = "{}/btk_summary_table_full.tsv".format(output_folder)
    counts_df_output_path = "{}/{}_phylum_counts_and_coverage.csv".format(
        output_folder, sample_id
    )

    if os.path.isdir(output_folder) == False:
        sys.stderr.write(
            "The directory with the output tables of the pipeline ({}) was not found\n".format(
                output_folder
            )
        )
        sys.exit(1)

    if os.path.isfile(main_results_table_path) == False:
        sys.stderr.write(
            "The main results table file ({}) was not found\n".format(
                main_results_table_path
            )
        )
        sys.exit(1)

    # Continue even if the BlobToolKit dataset is not found
    if os.path.isfile(btk_results_table_path) == False:
        sys.stderr.write(
            "Warning: BlobToolKit dataset ({}) was not found. Continuing without it.\n".format(
                btk_results_table_path
            )
        )

    main_df = None
    main_df = pd.read_csv(main_results_table_path)
    if main_df.shape[0] == 0:
        sys.stderr.write(
            "No rows were found in cobiont check results table ({})\n".format(
                main_results_table_path
            )
        )
        sys.exit(1)

    # Check for taxonomy columns in order of preference
    taxonomy_column = None
    if "buscoregions_phylum" in main_df.columns:
        taxonomy_column = "buscoregions_phylum"
    elif "buscogenes_phylum" in main_df.columns:
        taxonomy_column = "buscogenes_phylum"
    elif "btk_bestsum_phylum" in main_df.columns:
        taxonomy_column = "btk_bestsum_phylum"
    elif "fcs_gx_phylum" in main_df.columns:
        taxonomy_column = "fcs_gx_phylum"
    elif "nt_kraken_phylum" in main_df.columns:
        taxonomy_column = "nt_kraken_phylum"
    elif "tiara_classif" in main_df.columns:
        taxonomy_column = "tiara_classif"
    
    if taxonomy_column is None:
        sys.stderr.write(
            "No taxonomy column (buscoregions_phylum, buscogenes_phylum, btk_bestsum_phylum, fcs_gx_phylum, nt_kraken_phylum, or tiara_classif) was found in results table ({})\n".format(
                main_results_table_path
            )
        )
        sys.exit(1)
    
    # Check if coverage column exists
    if "coverage" not in main_df.columns:
        sys.stderr.write(
            "Column 'coverage' was not found in results table ({})\n".format(
                main_results_table_path
            )
        )
        sys.exit(1)
    
    # Check if length column exists
    if "length" not in main_df.columns:
        # If length column doesn't exist, add NaN values to indicate missing data
        sys.stderr.write(
            "Column 'length' was not found in results table ({}). Using NaN to indicate missing values.\n".format(
                main_results_table_path
            )
        )
        main_df["length"] = np.nan

    df = main_df
    df["merged_classif"] = df[taxonomy_column]
    df["merged_classif_source"] = taxonomy_column

    if "nt_kraken_phylum" in df.columns:
        ind = list(np.where(df["merged_classif"] == "no-hit")[0])
        ind2 = list(np.where(df["nt_kraken_phylum"].isna())[0])
        ind3 = [n for n in ind if n not in ind2]

        df.iloc[ind3, df.columns.get_loc("merged_classif")] = df.iloc[
            ind3, df.columns.get_loc("nt_kraken_phylum")
        ]
        df.iloc[ind3, df.columns.get_loc("merged_classif_source")] = "nt_kraken_phylum"

    if "tiara_classif" in df.columns:
        tiara_ind = list(np.where(df["merged_classif"] == "no-hit")[0])
        tiara_ind2 = list(np.where(df["tiara_classif"].isna())[0])
        tiara_ind3 = list(np.where(df["tiara_classif"] == "unknown")[0])
        tiara_ind = [
            n for n in tiara_ind if n not in tiara_ind2 and n not in tiara_ind3
        ]
        df.iloc[tiara_ind, df.columns.get_loc("merged_classif")] = df.iloc[
            tiara_ind, df.columns.get_loc("tiara_classif")
        ]
        df.iloc[tiara_ind, df.columns.get_loc("merged_classif_source")] = "tiara"

        df["merged_classif"] = df["merged_classif"].replace(
            "bacteria", "Bacteria-undef"
        )
        df["merged_classif"] = df["merged_classif"].replace(
            "eukarya", "Eukaryota-undef"
        )
        df["merged_classif"] = df["merged_classif"].replace(
            "prokarya", "Prokaryota-undef"
        )
        df["merged_classif"] = df["merged_classif"].replace("archaea", "Archaea-undef")

    # Write the modified DataFrame directly back to the original file
    # This is simpler and avoids issues with file operations
    df.to_csv(main_results_table_path, index=False)
    
    # Generate the counts file separately
    generate_counts_df(df, counts_df_output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output_folder",
        type=str,
        help="Path to the directory with the output tables of the pipeline",
    )
    parser.add_argument("sample_id", type=str, help="Sample ID")
    args = parser.parse_args()
    main(args.output_folder, args.sample_id)
