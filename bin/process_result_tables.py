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
    Creates a table that shows the number of sequences, mean coverage and sequence span per phylum
    """
    merged_classif_counts = df["merged_classif"].value_counts(dropna=False)
    cov_list = list()
    span_list = list()
    for classif_item, _ in merged_classif_counts.iteritems():
        ind = list(np.where(df["merged_classif"] == classif_item)[0])
        cov_values = df.iloc[ind, df.columns.get_loc("coverage")]
        length_values = df.iloc[ind, df.columns.get_loc("length")]
        cov_list.append(round(np.mean(cov_values), 2))
        span_sum = sum(length_values) / 1000000
        span_list.append(round(span_sum, 2))

    counts_df = merged_classif_counts.to_frame()
    counts_df["mean_coverage"] = cov_list
    counts_df["span_mb"] = span_list

    counts_df.rename(columns={"merged_classif": "number_of_sequences"}, inplace=True)
    counts_df.index.name = "phylum"
    counts_df.to_csv(counts_df_output_path, index=True)


def main(output_folder, sample_id):
    main_results_table_path = "{}/{}_contamination_check_merged_table.csv".format(output_folder, sample_id)
    btk_results_table_path = "{}/btk_summary_table_full.tsv".format(output_folder)
    output_df_path = "{}/{}_contamination_check_merged_table_extended.csv".format(output_folder, sample_id)
    counts_df_output_path = "{}/{}_phylum_counts_and_coverage.csv".format(output_folder, sample_id)

    if os.path.isdir(output_folder) == False:
        sys.stderr.write(
            "The directory with the output tables of the pipeline ({}) was not found\n".format(output_folder)
        )
        sys.exit(1)

    if os.path.isfile(main_results_table_path) == False:
        sys.stderr.write("The main results table file ({}) was not found\n".format(main_results_table_path))
        sys.exit(1)

    if os.path.isfile(btk_results_table_path) == False:
        sys.stderr.write(
            "Skipping the exporting of extended results table because the BlobToolKit dataset ({}) was not found\n".format(
                btk_results_table_path
            )
        )
        sys.exit(0)

    main_df = None
    main_df = pd.read_csv(main_results_table_path)
    if main_df.shape[0] == 0:
        sys.stderr.write("No rows were found in cobiont check results table ({})\n".format(main_results_table_path))
        sys.exit(1)

    if "btk_bestsum_phylum" not in main_df.columns:
        sys.stderr.write(
            "Column 'btk_bestsum_phylum' was not found in results table ({})\n".format(main_results_table_path)
        )
        sys.exit(1)

    df = main_df
    df["merged_classif"] = df["btk_bestsum_phylum"]
    df["merged_classif_source"] = "btk_bestsum_phylum"

    if "nt_kraken_phylum" in df.columns:
        ind = list(np.where(df["merged_classif"] == "no-hit")[0])
        ind2 = list(np.where(df["nt_kraken_phylum"].isna())[0])
        ind3 = [n for n in ind if n not in ind2]

        df.iloc[ind3, df.columns.get_loc("merged_classif")] = df.iloc[ind3, df.columns.get_loc("nt_kraken_phylum")]
        df.iloc[ind3, df.columns.get_loc("merged_classif_source")] = "nt_kraken_phylum"

    if "tiara_classif" in df.columns:
        tiara_ind = list(np.where(df["merged_classif"] == "no-hit")[0])
        tiara_ind2 = list(np.where(df["tiara_classif"].isna())[0])
        tiara_ind3 = list(np.where(df["tiara_classif"] == "unknown")[0])
        tiara_ind = [n for n in tiara_ind if n not in tiara_ind2 and n not in tiara_ind3]
        df.iloc[tiara_ind, df.columns.get_loc("merged_classif")] = df.iloc[
            tiara_ind, df.columns.get_loc("tiara_classif")
        ]
        df.iloc[tiara_ind, df.columns.get_loc("merged_classif_source")] = "tiara"

        df["merged_classif"] = df["merged_classif"].replace("bacteria", "Bacteria-undef")
        df["merged_classif"] = df["merged_classif"].replace("eukarya", "Eukaryota-undef")
        df["merged_classif"] = df["merged_classif"].replace("prokarya", "Prokaryota-undef")
        df["merged_classif"] = df["merged_classif"].replace("archaea", "Archaea-undef")

    df.to_csv(output_df_path, index=False)
    os.remove(main_results_table_path)
    # os.rename(output_df_path, main_results_table_path) To remove Nextflow confussion
    generate_counts_df(df, counts_df_output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_folder", type=str, help="Path to the directory with the output tables of the pipeline")
    parser.add_argument("sample_id", type=str, help="ToL ID of the sample")
    args = parser.parse_args()
    main(args.output_folder, args.sample_id)
