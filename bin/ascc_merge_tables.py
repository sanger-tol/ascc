#!/usr/bin/env python3

VERSION = "2.0.0"
DESCRIPTION = """
Script for merging contaminant check results into one table
Version: {VERSION}
---
Written by Eerik Aunin

Re-Written by Damon-Lee Pointon (dp24/DLBPointon)
"""

import argparse
import pandas as pd
import textwrap
import os
import sys
import general_purpose_functions as gpf


def parse_args():
    parser = argparse.ArgumentParser(
        prog="AsccMergeTables",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(DESCRIPTION),
    )
    parser.add_argument("-gc", "--gc_cov", required=True, type=str, help="GC Coverage file")
    parser.add_argument("-c", "--coverage", type=str, help="Coverage file")
    parser.add_argument("-t", "--tiara", type=str, help="Tiara file")
    parser.add_argument("-bk", "--bacterial_kraken", type=str, help="Bacterial Kraken file")
    parser.add_argument("-nk", "--nt_kraken", type=str, help="NT Kraken file")
    parser.add_argument("-nb", "--nt_blast", type=str, help="NT Blast file")
    parser.add_argument("-dr", "--dim_reduction_embeddings", type=str, help="Dimensional Reduction file")
    parser.add_argument("-nd", "--nr_diamond", type=str, help="NR Diamond file")
    parser.add_argument("-ud", "--uniprot_diamond", type=str, help="Uniprot Diamond file")
    parser.add_argument("-cv", "--contigviz", type=str, help="Contigviz file")
    parser.add_argument("-btk", "--blobtoolkit", type=str, help="Blobtoolkit file")
    parser.add_argument("-bb", "--busco_btk", type=str, help="Busco Blobtoolkit file")
    parser.add_argument("-fg", "--fcs_gx", type=str, help="FCS_GX file")
    parser.add_argument("-n", "--sample_name", type=str, help="Name for the sample")
    parser.add_argument("-m", "--markerscan", type=str, help="MarkerScan file")
    parser.add_argument("-v", "--version", action="version", version=VERSION)
    return parser.parse_args()


def check_paths(paths_dict, required_files):
    """
    Checks if a required file exists and exits with an error message if it doesn't
    """
    out_dict = dict()

    for data_type, input_file in paths_dict.items():
        if input == None:
            pass
        else:
            out_dict[data_type] = input_file

    return out_dict


def load_and_merge_dataframes(paths_dict):
    """
    Loads the tables with individual variables (GC content, coverage, kmer counts etc) and combines them into one table
    """
    gc_path = paths_dict["gc_content"]
    df = pd.read_csv(gc_path, sep="\t", header=None)
    if df.shape[0] > 0:
        df.columns = ["scaff", "gc"]
        df["gc"] = df["gc"] * 100
    else:
        sys.stderr.write("No rows were found in the GC content table ({})\n".format(gc_path))
        sys.exit(1)

    coverage_df = None
    if paths_dict["coverage"] is not None:
        coverage_df = pd.read_csv(paths_dict["coverage"], sep=",", header=None)
        if coverage_df.shape[0] > 0:
            coverage_df.columns = ["scaff", "coverage"]
        else:
            sys.stderr.write(f"No rows were found in the coverages table ({paths_dict['coverage']})\n")
            coverage_df = None

    tiara_df = None
    if paths_dict["tiara"] is not None:
        tiara_df = pd.read_csv(paths_dict["tiara"], sep="\t")
        if tiara_df.shape[0] > 0:
            tiara_df["tiara_classif"] = tiara_df["class_fst_stage"]
            tiara_snd_stage_hits = tiara_df.index[tiara_df["class_snd_stage"].notnull()]
            tiara_df["tiara_classif"][tiara_snd_stage_hits] = tiara_df["class_snd_stage"][tiara_snd_stage_hits]
            tiara_df = tiara_df.iloc[:, [0, 3]]
            tiara_df.columns = ["scaff", "tiara_classif"]
        else:
            sys.stderr.write("No rows were found in Tiara output table ({})\n".format(paths_dict["tiara"]))
            tiara_df = None

    bacterial_kraken_df = None
    if paths_dict["bacterial_kraken"] is not None:
        bacterial_kraken_df = pd.read_csv(paths_dict["bacterial_kraken"], sep=",")
        if bacterial_kraken_df.shape[0] > 0:
            bacterial_kraken_df.rename(columns={bacterial_kraken_df.columns[0]: "scaff"}, inplace=True)
            bacterial_kraken_df.rename(columns={"taxid": "kraken_taxid"}, inplace=True)
        else:
            sys.stderr.write(
                "No rows were found in bacterial Kraken output table ({})\n".format(paths_dict["bacterial_kraken"])
            )
            bacterial_kraken_df = None

    nt_kraken_df = None
    if paths_dict["nt_kraken"] is not None:
        nt_kraken_df = pd.read_csv(paths_dict["nt_kraken"], sep=",")
        if nt_kraken_df.shape[0] > 0:
            nt_kraken_df.rename(columns={nt_kraken_df.columns[0]: "scaff"}, inplace=True)
            nt_kraken_df.rename(columns={"taxid": "kraken_taxid"}, inplace=True)
        else:
            sys.stderr.write("No rows were found in nt Kraken output table ({})\n".format(paths_dict["nt_kraken"]))
            nt_kraken_df = None

    dim_reduction_df = None
    if paths_dict["dim_reduction_embeddings"] is not None:
        dim_reduction_df = pd.read_csv(paths_dict["dim_reduction_embeddings"], sep=",")
        if dim_reduction_df.shape[0] == 0:
            sys.stderr.write(
                "No rows were found in kmers dimensionality reduction output table ({})\n".format(
                    paths_dict["dim_reduction_embeddings"]
                )
            )
            dim_reduction_df = None

    btk_df = None
    if paths_dict["blobtoolkit"] is not None:
        btk_df = pd.read_csv(paths_dict["blobtoolkit"], header=0, delimiter="\t")
        if btk_df.shape[0] == 0:
            sys.stderr.write(
                "No rows were found in the BlobToolKit results table ({})\n".format(paths_dict["blobtoolkit"])
            )
            sys.exit(1)
        btk_renaming_dict = {"identifiers": "scaff", "bestsum_phylum": "btk_bestsum_phylum"}
        if "mapped_hifi_reads_sorted_cov" in btk_df.columns:
            btk_renaming_dict["mapped_hifi_reads_sorted_cov"] = "btk_cov"
        if "bestsum_phylum" in btk_df.columns:
            btk_renaming_dict["bestsum_phylum"] = "btk_bestsum_phylum"
        # {"identifiers": "scaff", "mapped_hifi_reads_sorted_cov": "btk_cov", "bestsum_phylum": "btk_bestsum_phylum"}

        btk_df.rename(columns=btk_renaming_dict, inplace=True)

        btk_selected_cols = [
            col for col in btk_df.columns if col in ["scaff", "length", "btk_cov", "btk_bestsum_phylum"]
        ]
        if len(btk_selected_cols) > 0:
            btk_df = btk_df[btk_selected_cols]
        else:
            btk_df = None

    btk_busco_df = None
    if paths_dict["btk_busco"] is not None:
        btk_busco_df = pd.read_csv(paths_dict["btk_busco"], header=0, delimiter="\t")
        if btk_busco_df.shape[0] == 0:
            sys.stderr.write(
                "No rows were found in the BUSCO-based BlobToolKit results table ({})\n".format(paths_dict["btk_busco"])
            )
            sys.exit(1)
        btk_busco_renaming_dict = {"identifiers": "scaff"}

        btk_busco_df.rename(columns=btk_busco_renaming_dict, inplace=True)

        btk_busco_selected_cols = [
            col
            for col in btk_busco_df.columns
            if col
            in [
                "scaff",
                "buscogenes_superkingdom",
                "buscogenes_kingdom",
                "buscogenes_phylum",
                "buscogenes_class",
                "buscogenes_order",
                "buscogenes_family",
                "buscogenes_genus",
                "buscogenes_species",
                "buscoregions_superkingdom",
                "buscoregions_kingdom",
                "buscoregions_phylum",
                "buscoregions_class",
                "buscoregions_order",
                "buscoregions_family",
                "buscoregions_genus",
                "buscoregions_species",
            ]
        ]
        if len(btk_busco_selected_cols) > 0:
            btk_busco_df = btk_busco_df[btk_busco_selected_cols]
        else:
            btk_busco_df = None

    fcs_gx_df = None
    if paths_dict["fcs_gx"] is not None:
        fcs_gx_df = pd.read_csv(paths_dict["fcs_gx"], sep=",")
        if fcs_gx_df.shape[0] == 0:
            sys.stderr.write("No rows were found in FCS-GX output table ({})\n".format(paths_dict["fcs_gx"]))
            fcs_gx_df = None

    nt_blast_df = None
    if paths_dict["nt_blast"] is not None:
        nt_blast_df = pd.read_csv(paths_dict["nt_blast"], sep=",")
        if nt_blast_df.shape[0] == 0:
            sys.stderr.write("No rows were found in nt BLAST output table ({})\n".format(paths_dict["nt_blast"]))
            nt_blast_df = None

    nr_diamond_df = None
    if paths_dict["nr_diamond"] is not None:
        nr_diamond_df = pd.read_csv(paths_dict["nr_diamond"], sep=",")
        if nr_diamond_df.shape[0] == 0:
            sys.stderr.write("No rows were found in nr Diamond output table ({})\n".format(paths_dict["nr_diamond"]))
            nr_diamond_df = None

    uniprot_diamond_df = None
    if paths_dict["uniprot_diamond"] is not None:
        uniprot_diamond_df = pd.read_csv(paths_dict["uniprot_diamond"], sep=",")
        if uniprot_diamond_df.shape[0] == 0:
            sys.stderr.write(
                "No rows were found in Uniprot Diamond output table ({})\n".format(paths_dict["uniprot_diamond"])
            )
            uniprot_diamond_df = None

    cobiontid_markerscan_df = None
    if paths_dict["cobiontid_markerscan"] is not None:
        cobiontid_markerscan_df = pd.read_csv(paths_dict["cobiontid_markerscan"], sep=",")
        if cobiontid_markerscan_df.shape[0] == 0:
            sys.stderr.write(
                "No rows were found in CobiontID MarkerScan output table ({})\n".format(
                    paths_dict["cobiontid_markerscan"]
                )
            )
            uniprot_diamond_df = None

    contigviz_df = None
    if paths_dict["contigviz"] is not None:
        contigviz_df = pd.read_csv(paths_dict["contigviz"], sep=",")
        if contigviz_df.shape[0] == 0:
            sys.stderr.write("No rows were found in ContigViz output table ({})\n".format(paths_dict["contigviz"]))
            contigviz_df = None

    if coverage_df is not None:
        df = pd.merge(df, coverage_df, on="scaff", how="outer")
    if tiara_df is not None:
        df = pd.merge(df, tiara_df, on="scaff", how="outer")
    if bacterial_kraken_df is not None:
        df = pd.merge(df, bacterial_kraken_df, on="scaff", how="outer")
    if nt_kraken_df is not None:
        df = pd.merge(df, nt_kraken_df, on="scaff", how="outer")
    if dim_reduction_df is not None:
        df = pd.merge(df, dim_reduction_df, on="scaff", how="outer")
    if nt_blast_df is not None:
        df = pd.merge(df, nt_blast_df, on="scaff", how="outer")
    if nr_diamond_df is not None:
        df = pd.merge(df, nr_diamond_df, on="scaff", how="outer")
    if uniprot_diamond_df is not None:
        df = pd.merge(df, uniprot_diamond_df, on="scaff", how="outer")
    if fcs_gx_df is not None:
        df = pd.merge(df, fcs_gx_df, on="scaff", how="outer")
    if cobiontid_markerscan_df is not None:
        df = pd.merge(df, cobiontid_markerscan_df, on="scaff", how="outer")
    if contigviz_df is not None:
        df = pd.merge(df, contigviz_df, on="scaff", how="outer")
    if btk_df is not None:
        df = pd.merge(df, btk_df, on="scaff", how="outer")
    if btk_busco_df is not None:
        df = pd.merge(df, btk_busco_df, on="scaff", how="outer")

    return df


def main(args):
    paths_dict = dict()
    paths_dict["gc_content"] = args.gc_cov
    paths_dict["coverage"] = args.coverage
    paths_dict["tiara"] = args.tiara
    paths_dict["bacterial_kraken"] = args.bacterial_kraken
    paths_dict["nt_kraken"] = args.nt_kraken
    paths_dict["nt_blast"] = args.nt_blast
    paths_dict["dim_reduction_embeddings"] = args.dim_reduction_embeddings
    paths_dict["nr_diamond"] = args.nr_diamond
    paths_dict["uniprot_diamond"] = args.uniprot_diamond
    paths_dict["cobiontid_markerscan"] = args.markerscan
    paths_dict["contigviz"] = args.contigviz
    paths_dict["blobtoolkit"] = args.blobtoolkit
    paths_dict["btk_busco"] = args.busco_btk
    paths_dict["fcs_gx"] = args.fcs_gx

    required_files = ["gc_content"]

    paths_dict = check_paths(paths_dict, required_files)
    df = load_and_merge_dataframes(paths_dict)
    df.to_csv(f"{args.sample_name}_contamination_check_merged_table.csv", index=False)

    if (
        paths_dict["nt_blast"]
        and paths_dict["nr_diamond"]
        and paths_dict["uniprot_diamond"]
        and paths_dict["coverage"]
        and paths_dict["tiara"]
        and paths_dict["nt_kraken"]
    ):
        process_results_tables_command = f"process_result_tables.py . {args.sample_name}"
        gpf.run_system_command(process_results_tables_command)
    else:
        sys.stderr.write(
            f"Skipping generating the {args.sample_name}_phylum_counts_and_coverage.csv file, as the variables used in this run do not include all the required variables for this (nt_blast, nr_diamond, uniprot_diamond, coverage, tiara, nt_kraken)\n"
        )


if __name__ == "__main__":
    main(parse_args())
