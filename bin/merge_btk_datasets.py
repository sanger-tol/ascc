#!/usr/bin/env python3

VERSION = "2.0.0"
DESCRIPTION = f"""
---
Script for merging BlobToolKit datasets from the createBTKdatasets output directory.
Version: {VERSION}
---

Written by Eerik Aunin (ea10)

Modified by Damon-Lee Pointon (@dp24/@DLBPointon)

"""

import json
from pathlib import Path
import shutil
import os
import sys
import argparse
import textwrap
import general_purpose_functions as gpf
import gzip


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="mergeBTKdatasets",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(DESCRIPTION),
    )
    parser.add_argument(
        "-m", "--main_btk_datasets", required=True, type=str, help="The btk_datasets generated by createBTKdatasets"
    )
    parser.add_argument(
        "-b",
        "--btk_busco_datasets",
        type=str,
        help="Path to the BTK dataset (blobdir) created by the BUSCO-based BTK pipeline",
    )
    parser.add_argument(
        "-s",
        "--btk_busco_summary_full",
        type=str,
        help="The btk_datasets generated by createBTKdatasets",
    )
    parser.add_argument(
        "-o",
        "--new_output_directory",
        default="merged_datasets",
        type=str,
        help="The new output directory for the merged datasets",
    )
    parser.add_argument("-v", "--version", action="version", version=VERSION)

    return parser.parse_args(argv)


def load_json(filename):
    """
    Loads a JSON file and returns it as a dictionary
    """
    json_contents = None
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt", encoding="UTF-8") as zipfile:
            json_contents = json.load(zipfile)
    else:
        with open(filename) as f:
            json_contents = json.load(f)
    return json_contents


def create_meta_json_contents(main_btk_dataset_folder, btk_busco_dataset_folder):
    """
    Creates the contents for the meta.json file for the new BTK dataset by combining the two meta.json files from the input directories
    """
    for folder in (main_btk_dataset_folder, btk_busco_dataset_folder):
        if os.path.isdir(folder) is False:
            sys.stderr.write(
                f"Skipping the merging of the main BTK dataset and the BUSCO-based BTK dataset, as directory {folder} was not found)\n"
            )
            sys.exit(0)

    main_btk_json_path = f"{main_btk_dataset_folder}/meta.json"
    btk_busco_json_path = f"{btk_busco_dataset_folder}/meta.json.gz"
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
    return merged_dict


def detect_buscogenes_variables(merged_jsons_dict):
    """
    Goes through the content of merged meta.json file (derived from both BTK datasets) and detects if buscogenes
        variables are present
    """
    buscogenes_present_flag = False
    fields = merged_jsons_dict["fields"]
    for field in fields:
        field_id = field["id"]
        if field_id == "taxonomy":
            for item in field["children"]:
                if item["id"] == "buscogenes":
                    buscogenes_present_flag = True
                    break
    return buscogenes_present_flag


def main(args):
    if os.path.isdir(args.main_btk_datasets) is False:
        sys.stderr.write(f"The BlobToolKit dataset ({args.main_btk_datasets}) was not found!\n")
        sys.exit(1)

    if os.path.isdir(args.btk_busco_datasets) is False:
        sys.stderr.write(
            f"The blobdir of BUSCO-based BlobToolKit Snakemake pipeline run does not exist at {args.btk_busco_datasets}, skipping the merging of BTK datasets\n"
        )
        sys.exit(0)

    not_copying_list = [
        "identifiers.json.gz",
        "gc_data.json.gz",
        "length_data.json.gz",
        "ncount_data.json.gz",
        "meta.json.gz",
    ]

    Path(args.new_output_directory).mkdir(parents=True, exist_ok=True)

    main_btk_dataset_files = [
        f for f in os.listdir(args.main_btk_datasets) if os.path.isfile(os.path.join(args.main_btk_datasets, f))
    ]
    main_btk_dataset_files = [f for f in main_btk_dataset_files if f not in not_copying_list]
    for main_btk_dataset_file in main_btk_dataset_files:
        main_btk_dataset_file_full_path = f"{args.main_btk_datasets}/{main_btk_dataset_file}"
        copied_file_full_path = os.path.abspath(f"{args.new_output_directory}/{main_btk_dataset_file}")
        shutil.copy(main_btk_dataset_file_full_path, copied_file_full_path)

    btk_busco_files = [
        f for f in os.listdir(args.btk_busco_datasets) if os.path.isfile(os.path.join(args.btk_busco_datasets, f))
    ]
    for btk_busco_file in btk_busco_files:
        btk_busco_file_full_path = f"{args.btk_busco_datasets}/{btk_busco_file}"
        copied_file_full_path = os.path.abspath(f"{args.new_output_directory}/{btk_busco_file}")
        shutil.copy(btk_busco_file_full_path, copied_file_full_path)

    merged_jsons_dict = create_meta_json_contents(args.main_btk_datasets, args.btk_busco_datasets)
    meta_json_outpath = f"{args.new_output_directory}/meta.json"

    with open(meta_json_outpath, "w") as json_outfile:
        json.dump(merged_jsons_dict, json_outfile, indent=1, sort_keys=True)

    buscogenes_present_flag = detect_buscogenes_variables(merged_jsons_dict)

    btk_busco_table_outpath = f"{args.new_output_directory}/btk_busco_summary_table_full.tsv"

    btk_busco_table_exporting_command = f"blobtools filter --table {btk_busco_table_outpath} --table-fields identifiers,buscoregions_superkingdom,buscoregions_kingdom,buscoregions_phylum,buscoregions_class,buscoregions_order,buscoregions_family,buscoregions_genus,buscoregions_species"
    if buscogenes_present_flag == True:
        btk_busco_table_exporting_command += ",buscogenes_superkingdom,buscogenes_kingdom,buscogenes_phylum,buscogenes_class,buscogenes_order,buscogenes_family,buscogenes_genus,buscogenes_species"
    btk_busco_table_exporting_command += f" {args.new_output_directory}"

    gpf.run_system_command(btk_busco_table_exporting_command)


if __name__ == "__main__":
    main(parse_args())
