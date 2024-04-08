#!/usr/bin/env python3
"""
Script for merging contaminant check results into one table
"""

import argparse
import pandas as pd
import os
import sys
import general_purpose_functions as gpf


def check_paths(paths_dict, required_files):
    """
    Checks if a required file exists and exits with an error message if it doesn't
    """
    out_dict = dict()
    for data_type, input_file_path in paths_dict.items():
        out_dict[data_type] = None
        if os.path.isfile(input_file_path) == False:
            if data_type in required_files:
                sys.stderr.write("Input file {} was not found\n".format(input_file_path))
                sys.exit(1)
        else:
            if os.stat(input_file_path).st_size == 0:
                sys.stderr.write("Warning: the file {} is empty and will therefore not be included in the final results\n".format(input_file_path))
            else:
                out_dict[data_type] = input_file_path
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
    coverage_file_path = paths_dict["coverage"]
    if coverage_file_path is not None:
        if os.stat(coverage_file_path).st_size > 0:
            coverage_df = pd.read_csv(coverage_file_path, sep=",", header=None)
            if coverage_df.shape[0] > 0:
                coverage_df.columns = ["scaff", "coverage"]
            else:
                sys.stderr.write("No rows were found in the coverages table ({})\n".format(coverage_file_path))
                coverage_df = None
        else:
            sys.stderr.write("Warning: the output file for PacBio coverage ({}) is empty\n".format(coverage_file_path))

    tiara_df = None
    if paths_dict["tiara"] is not None:
        tiara_df = pd.read_csv(paths_dict["tiara"], sep="\t")
        if tiara_df.shape[0] > 0:
            tiara_df["tiara_classif"] = tiara_df["class_fst_stage"]
            tiara_snd_stage_hits = tiara_df.index[tiara_df["class_snd_stage"].notnull()]
            tiara_df["tiara_classif"][tiara_snd_stage_hits] = tiara_df["class_snd_stage"][tiara_snd_stage_hits]
            tiara_df = tiara_df.iloc[:,[0, 3]]
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
            sys.stderr.write("No rows were found in bacterial Kraken output table ({})\n".format(paths_dict["bacterial_kraken"]))
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
            sys.stderr.write("No rows were found in kmers dimensionality reduction output table ({})\n".format(paths_dict["dim_reduction_embeddings"]))
            dim_reduction_df = None

    btk_df = None

    if paths_dict["blobtoolkit"] is not None:
        btk_df = pd.read_csv(paths_dict["blobtoolkit"], header=0, delimiter="\t")
        if btk_df.shape[0] == 0:
            sys.stderr.write("No rows were found in the BlobToolKit results table ({})\n".format(paths_dict["blobtoolkit"]))
            sys.exit(1)
        btk_renaming_dict = {"identifiers": "scaff", "bestsum_phylum": "btk_bestsum_phylum"}
        if "mapped_hifi_reads_sorted_cov" in btk_df.columns:
            btk_renaming_dict["mapped_hifi_reads_sorted_cov"] = "btk_cov"
        if "bestsum_phylum" in btk_df.columns:
            btk_renaming_dict["bestsum_phylum"] = "btk_bestsum_phylum"
        #{"identifiers": "scaff", "mapped_hifi_reads_sorted_cov": "btk_cov", "bestsum_phylum": "btk_bestsum_phylum"}

        btk_df.rename(columns = btk_renaming_dict, inplace=True)

        btk_selected_cols = [col for col in btk_df.columns if col in ["scaff", "length", "btk_cov", "btk_bestsum_phylum"]]
        if len(btk_selected_cols) > 0:
            btk_df = btk_df[btk_selected_cols]
        else:
            btk_df = None

    btk_busco_df = None

    if paths_dict["btk_busco"] is not None:
        btk_busco_df = pd.read_csv(paths_dict["btk_busco"], header=0, delimiter="\t")
        if btk_busco_df.shape[0] == 0:
            sys.stderr.write("No rows were found in the BUSCO-based BlobToolKit results table ({})\n".format(paths_dict["btk_busco"]))
            sys.exit(1)
        btk_busco_renaming_dict = {"identifiers": "scaff"}
        #if "mapped_hifi_reads_sorted_cov" in btk_df.columns:
        #    btk_renaming_dict["mapped_hifi_reads_sorted_cov"] = "btk_cov"
        #if "bestsum_phylum" in btk_df.columns:
        #    btk_renaming_dict["bestsum_phylum"] = "btk_bestsum_phylum"
        #{"identifiers": "scaff", "mapped_hifi_reads_sorted_cov": "btk_cov", "bestsum_phylum": "btk_bestsum_phylum"}

        btk_busco_df.rename(columns = btk_busco_renaming_dict, inplace=True)

        btk_busco_selected_cols = [col for col in btk_busco_df.columns if col in ["scaff", "buscogenes_superkingdom", "buscogenes_kingdom", "buscogenes_phylum", "buscogenes_class", "buscogenes_order", "buscogenes_family", "buscogenes_genus", "buscogenes_species", "buscoregions_superkingdom", "buscoregions_kingdom", "buscoregions_phylum", "buscoregions_class", "buscoregions_order", "buscoregions_family", "buscoregions_genus", "buscoregions_species"]]
        if len(btk_busco_selected_cols) > 0:
            btk_busco_df = btk_busco_df[btk_busco_selected_cols]
        else:
            btk_busco_df = None

        #df = pd.merge(main_df, btk_df, on="scaff", how="outer")


    #if paths_dict["blobtoolkit"] is not None:
    #    #if 'A' in df.columns:
    #    blobtoolkit_df = pd.read_csv(paths_dict["blobtoolkit"], header=0, delimiter="\t")
    #    if blobtoolkit_df.shape[0] > 0:
    #        blobtoolkit_df = blobtoolkit_df[["identifiers", "bestsum_phylum"]]
    #        blobtoolkit_df.columns = ["scaff", "btk_bestsum_phylum"]
    #    else:
    #        sys.stderr.write("No rows were found in BlobToolKit output table ({})\n".format(paths_dict["blobtoolkit"]))
    #        blobtoolkit_df = None

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
            sys.stderr.write("No rows were found in Uniprot Diamond output table ({})\n".format(paths_dict["uniprot_diamond"]))
            uniprot_diamond_df = None
    cobiontid_markerscan_df = None
    if paths_dict["cobiontid_markerscan"] is not None:
        cobiontid_markerscan_df = pd.read_csv(paths_dict["cobiontid_markerscan"], sep=",")
        if cobiontid_markerscan_df.shape[0] == 0:
            sys.stderr.write("No rows were found in CobiontID MarkerScan output table ({})\n".format(paths_dict["cobiontid_markerscan"]))
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


def main(data_folder, out_path, sample_name):
    paths_dict = dict()
    paths_dict["gc_content"] = "{}/gc.txt".format(data_folder)
    paths_dict["coverage"] = "{}/pacbio_reads_coverage.txt".format(data_folder)
    paths_dict["tiara"] = "{}/tiara_out.txt".format(data_folder)
    paths_dict["bacterial_kraken"] = "{}/bacterial_kraken_lineage.txt".format(data_folder)
    paths_dict["nt_kraken"] = "{}/nt_kraken_lineage.txt".format(data_folder)
    paths_dict["nt_blast"] = "{}/BLAST_results_with_lineage.csv".format(data_folder)
    paths_dict["dim_reduction_embeddings"] = "{}/kmers_dim_reduction_embeddings.csv".format(data_folder)
    paths_dict["nr_diamond"] = "{}/nr_diamond_blastx_top_hits.csv".format(data_folder)
    paths_dict["uniprot_diamond"] = "{}/uniprot_diamond_blastx_top_hits.csv".format(data_folder)
    paths_dict["cobiontid_markerscan"] = "{}/cobiontid_markerscan.csv".format(data_folder)
    paths_dict["contigviz"] = "{}/contigviz_results.csv".format(data_folder)
    paths_dict["blobtoolkit"] = "{}/btk_summary_table_full.tsv".format(data_folder)
    paths_dict["btk_busco"] = "{}/btk_busco_summary_table_full.tsv".format(data_folder)
    paths_dict["fcs_gx"] = "{}/fcs-gx_summary.csv".format(data_folder)

    required_files = ["gc_content"]

    paths_dict = check_paths(paths_dict, required_files)
    df = load_and_merge_dataframes(paths_dict)
    df.to_csv(out_path, index=False)

    if paths_dict["nt_blast"] is not None and paths_dict["nr_diamond"] is not None and paths_dict["uniprot_diamond"] is not None and paths_dict["coverage"] is not None and paths_dict["tiara"] is not None and paths_dict["nt_kraken"] is not None:
        process_results_tables_command = "process_result_tables.py {} {}".format(data_folder, sample_name)
        gpf.run_system_command(process_results_tables_command)
    else:
        sys.stderr.write("Skipping generating the {}_phylum_counts_and_coverage.csv file, as the variables used in this run do not include all the required variables for this (nt_blast, nr_diamond, uniprot_diamond, coverage, tiara, nt_kraken)\n".format(sample_name))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_folder", type=str, help="Path to folder with ASG contamination check result files of individual steps")
    parser.add_argument("out_path", type=str, help="Path for output CSV file")
    parser.add_argument("--sample_name", type=str, help="Sample name (e.g. ToLID)", default="unnamed")
    args = parser.parse_args()
    main(args.data_folder, args.out_path, args.sample_name)