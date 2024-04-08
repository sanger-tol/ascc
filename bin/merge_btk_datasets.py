#!/usr/bin/env python3
"""
Script for merging BTK datasets from the this pipeline and the BUSCO-based Snakemake BTK pipeline
"""

import json
from pathlib import Path
import shutil
import os
import sys
import argparse
import general_purpose_functions as gpf


def load_json(filename):
    """ Loads a JSON file and returns it as a dictionary """
    with open(filename) as f:
        return json.load(f)


def create_meta_json(main_btk_dataset_folder, btk_busco_dataset_folder, combined_dataset_folder):
    """
    Creates a meta.json file for the new BTK dataset by combining the two meta.json files from the input directories
    """
    for folder in (main_btk_dataset_folder, btk_busco_dataset_folder):
        if os.path.isdir(folder) is False:
            sys.stderr.write(f"Skipping the merging of the main BTK dataset and the BUSCO-based BTK dataset, as directory {folder} was not found)\n")
            sys.exit(0)

    main_btk_json_path = f"{main_btk_dataset_folder}/meta.json"
    btk_busco_json_path = f"{btk_busco_dataset_folder}/meta.json"
    for json_path in (main_btk_json_path, btk_busco_json_path):
        if os.path.isfile(json_path) is False:
            sys.stderr.write(f"File {json_path} not found)\n")
            sys.exit(1)
    main_meta_dict = load_json(main_btk_json_path)
    btk_busco_meta_dict = load_json(btk_busco_json_path)

    merged_dict = btk_busco_meta_dict.copy()

    keys_to_skip = []
    fields = main_meta_dict["fields"]
    for field in fields:
        field_id = field["id"]

        if field_id == "taxonomy":
            btk_main_taxonomy_field = field.copy()
            btk_main_taxonomy_field["id"] = "btk_main_taxonomy"
            btk_main_taxonomy_field["name"] = "btk_main_taxonomy"
            merged_dict["fields"].append(btk_main_taxonomy_field)
        else:
            if field_id not in keys_to_skip:
                merged_dict["fields"].append(field)


    meta_json_outpath = f"{combined_dataset_folder}/meta.json"
    with open(meta_json_outpath, "w") as json_outfile:
        json.dump(merged_dict, json_outfile, indent=1, sort_keys=True)


def main(main_btk_dataset_folder, btk_busco_dataset_folder, combined_dataset_folder, pipeline_output_folder, skip_renaming_folders):
    if os.path.isdir(main_btk_dataset_folder) is False:
        sys.stderr.write(f"The BlobToolKit dataset ({main_btk_dataset_folder}) was not found\n")
        sys.exit(1)

    if os.path.isdir(btk_busco_dataset_folder) is False:
        sys.stderr.write(f"The blobdir of BUSCO-based BlobToolKit Snakemake pipeline run does not exist at {btk_busco_dataset_folder}, skipping the merging of BTK datasets\n")
        sys.exit(0)

    not_copying_list = ["identifiers.json", "gc_data.json", "length_data.json", "ncount_data.json", "meta.json"]

    Path(combined_dataset_folder).mkdir(parents=True, exist_ok=True)

    main_btk_dataset_files = [f for f in os.listdir(main_btk_dataset_folder) if os.path.isfile(os.path.join(main_btk_dataset_folder, f))]
    main_btk_dataset_files = [f for f in main_btk_dataset_files if f not in not_copying_list]
    for main_btk_dataset_file in main_btk_dataset_files:
        main_btk_dataset_file_full_path = f"{main_btk_dataset_folder}/{main_btk_dataset_file}"
        copied_file_full_path = f"{combined_dataset_folder}/{main_btk_dataset_file}"
        shutil.copy(main_btk_dataset_file_full_path, copied_file_full_path)

    btk_busco_files = [f for f in os.listdir(btk_busco_dataset_folder) if os.path.isfile(os.path.join(btk_busco_dataset_folder, f))]
    for btk_busco_file in btk_busco_files:
        btk_busco_file_full_path = f"{btk_busco_dataset_folder}/{btk_busco_file}"
        copied_file_full_path = f"{combined_dataset_folder}/{btk_busco_file}"
        shutil.copy(btk_busco_file_full_path, copied_file_full_path)

    create_meta_json(main_btk_dataset_folder, btk_busco_dataset_folder, combined_dataset_folder)
    old_main_btk_dataset_folder = main_btk_dataset_folder + "_without_busco"

    if skip_renaming_folders is False:
        os.rename(main_btk_dataset_folder, old_main_btk_dataset_folder)
        os.rename(combined_dataset_folder, main_btk_dataset_folder)

    btk_busco_table_outpath = f"{pipeline_output_folder}/btk_busco_summary_table_full.tsv"
    btk_busco_table_exporting_command = f"blobtools filter --table {btk_busco_table_outpath} --table-fields identifiers,buscogenes_superkingdom,buscogenes_kingdom,buscogenes_phylum,buscogenes_class,buscogenes_order,buscogenes_family,buscogenes_genus,buscogenes_species,buscoregions_superkingdom,buscoregions_kingdom,buscoregions_phylum,buscoregions_class,buscoregions_order,buscoregions_family,buscoregions_genus,buscoregions_species {main_btk_dataset_folder}"
    gpf.run_system_command(btk_busco_table_exporting_command)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("main_btk_dataset_folder", type=str, help="Path to the BTK dataset (blobdir) created from the output of the steps of this pipeline")
    parser.add_argument("btk_busco_dataset_folder", type=str, help="Path to the BTK dataset (blobdir) created by the BUSCO-based Snakemake BTK pipeline")
    parser.add_argument("combined_dataset_folder", type=str, help="Path for creating a new BTK dataset (blobdir) that combines the two input BTK datasets")
    parser.add_argument("pipeline_output_folder", type=str, help="Path to the directory with the output tables of the pipeline")
    parser.add_argument("--skip_renaming_folders", dest="skip_renaming_folders", help="Optional boolean argument. If set to true, the script skips the renaming of the input BTK dataset directories after creating the merged BTK dataset", action="store_true")
    args = parser.parse_args()
    main(args.main_btk_dataset_folder, args.btk_busco_dataset_folder, args.combined_dataset_folder, args.pipeline_output_folder, args.skip_renaming_folders)