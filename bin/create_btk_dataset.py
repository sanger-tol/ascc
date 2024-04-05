#!/usr/bin/env python3
"""
Script for creating a BlobToolKit dataset
"""

import general_purpose_functions as gpf
import argparse
from pathlib import Path
import sys
import os.path


def create_assembly_yaml(assembly_yaml_path, assembly_alias, taxon_name):
    """
    Creates the assembly YAML file for creating a BlobToolKit dataset
    """
    if ".gz" in assembly_alias:
        assembly_alias = assembly_alias.replace(".gz", "_gz")
    out_string = "assembly:\n  accession: NA\n  alias: {}\n  record_type: scaffold\n  bioproject: NA\n  biosample: NA\ntaxon:\n  name: {}".format(assembly_alias, taxon_name)
    with open(assembly_yaml_path, "w") as f:
        f.write(out_string)


def tiara_results_to_btk_format(tiara_results_path, outfile_path):
    """
    Reformatting Tiara output file so that the summarised results of the first and second pass of Tiara can be
        added to a BlobToolKit dataset
    """
    tiara_data = gpf.l(tiara_results_path)
    tiara_data = tiara_data[1:len(tiara_data)]
    with open(outfile_path, "w") as f:
        f.write("identifier\ttiara\n")
        for line in tiara_data:
            split_line = line.split()
            if len(split_line) != 3:
                sys.stderr.write("Failed to parse the Tiara results file {}\n".format(tiara_results_path))
                sys.exit(1)
            first_pass_result = split_line[1]
            second_pass_result = split_line[2]
            if second_pass_result != "n/a":
                first_pass_result = second_pass_result
            f.write(split_line[0] + "\t" + first_pass_result + "\n")


def detect_dim_reduction_methods(kmers_dim_reduction_output_path):
    """
    Parses the header of the kmers dimensionality reduction report file to detect which dimensionality reduction methods were used
    """
    header_string = None
    with open(kmers_dim_reduction_output_path) as f:
        header_string = f.readline()
    header_string = header_string.strip()
    split_header = header_string.split(",")
    dim_reduction_methods = list()
    for header_item in split_header:
        if header_item.startswith("embedding_"):
            if header_item.startswith("embedding_x_"):
                header_item = header_item.split("embedding_x_")[1]
            elif header_item.startswith("embedding_y_"):
                header_item = header_item.split("embedding_y_")[1]
            if header_item not in dim_reduction_methods:
                dim_reduction_methods.append(header_item)
    return dim_reduction_methods


def add_custom_variables_to_btk_dataset(pipeline_run_folder, btk_dataset_folder):
    """
    Script for adding custom variables (e.g. Tiara results and PCA results) to the BlobToolKit dataset
    """
    pipeline_output_folder = pipeline_run_folder + "/collected_tables"
    if os.path.isdir(pipeline_output_folder) == False:
        sys.stderr.write("The directory for the output tables of the pipeline ({}) was not found\n".format(pipeline_output_folder))
        sys.exit(1)
    if os.path.isdir(btk_dataset_folder) == False:
        sys.stderr.write("The BlobToolKit dataset directory ({}) was not found\n".format(btk_dataset_folder))
        sys.exit(1)
    tiara_raw_output_path = pipeline_output_folder + "/tiara_out.txt"
    if os.path.isfile(tiara_raw_output_path) and os.stat(tiara_raw_output_path).st_size > 0:
        tiara_reformatted_output_path = pipeline_output_folder + "/tiara_out_btk_format.tsv"
        tiara_results_to_btk_format(tiara_raw_output_path, tiara_reformatted_output_path)
        add_tiara_command = 'blobtools add --text {} --text-delimiter "\t" --text-cols "identifier=identifiers,tiara=tiara" --text-header {}'.format(tiara_reformatted_output_path, btk_dataset_folder)
        gpf.run_system_command(add_tiara_command)

    kmers_dim_reduction_output_path = pipeline_output_folder + "/kmers_dim_reduction_embeddings.csv"
    if os.path.isfile(kmers_dim_reduction_output_path) and os.stat(kmers_dim_reduction_output_path).st_size > 0:
        used_dim_reduction_methods = detect_dim_reduction_methods(kmers_dim_reduction_output_path)
        for dim_reduction_method in used_dim_reduction_methods:
            add_embedding_command = 'blobtools add --text {path} --text-delimiter "," --text-cols scaff=identifiers,embedding_x_{dim_reduction_method}=embedding_x_{dim_reduction_method},embedding_y_{dim_reduction_method}=embedding_y_{dim_reduction_method} --text-header {btk_dataset_folder}'.format(path=kmers_dim_reduction_output_path, dim_reduction_method=dim_reduction_method, btk_dataset_folder=btk_dataset_folder)
            gpf.run_system_command(add_embedding_command)


    kraken_lineage_path = pipeline_output_folder + "/nt_kraken_lineage.txt"
    if os.path.isfile(kraken_lineage_path) and os.stat(kraken_lineage_path).st_size > 0:
        for taxonomy_level in ("species", "genus", "family", "order", "class", "phylum", "kingdom", "domain"):
            add_kraken_command = 'blobtools add --text {} --text-delimiter "," --text-cols scaff=identifiers,nt_kraken_{}=nt_kraken_{} --text-header {}'.format(kraken_lineage_path, taxonomy_level, taxonomy_level, btk_dataset_folder)
            gpf.run_system_command(add_kraken_command)

    fcs_gx_output_path = pipeline_output_folder + "/fcs-gx_summary.csv"
    if os.path.isfile(fcs_gx_output_path) and os.stat(fcs_gx_output_path).st_size > 0:
        add_fcs_gx_results_command = 'blobtools add --text {} --text-delimiter "," --text-cols "scaff=identifiers,fcs_gx_top_tax_name=fcs_gx_top_tax_name,fcs_gx_div=fcs_gx_div,fcs_gx_action=fcs_gx_action" --text-header {}'.format(fcs_gx_output_path, btk_dataset_folder)
        gpf.run_system_command(add_fcs_gx_results_command)

    #cobiontid_markerscan_json_file_path = run_folder + "/" + sample_id + ".json"
    #cobiontid_scaffs_json_to_csv(json_file_path, out_folder + "/cobiontid_markerscan.csv")
    cobiontid_markerscan_output_path = pipeline_output_folder + "/cobiontid_markerscan.csv"
    if os.path.isfile(cobiontid_markerscan_output_path) and os.stat(cobiontid_markerscan_output_path).st_size > 0:
        add_cobiontid_markerscan_results_command = 'blobtools add --text {} --text-delimiter "," --text-cols "scaff=identifiers,CobiontID_MarkerScan_embl_ebi_ena=CobiontID_MarkerScan_embl_ebi_ena,CobiontID_MarkerScan_slv=CobiontID_MarkerScan_slv,CobiontID_MarkerScan_Cluster=CobiontID_MarkerScan_Cluster" --text-header {}'.format(cobiontid_markerscan_output_path, btk_dataset_folder)
        gpf.run_system_command(add_cobiontid_markerscan_results_command)

    cobiontid_contigviz_output_path = pipeline_output_folder + "/contigviz_results.csv"
    if os.path.isfile(cobiontid_contigviz_output_path) and os.stat(cobiontid_contigviz_output_path).st_size > 0:
        add_cobiontid_contigviz_results_command = 'blobtools add --text {} --text-delimiter "," --text-cols "scaff=identifiers,ContigViz_UMAP1=ContigViz_UMAP1,ContigViz_UMAP2=ContigViz_UMAP2,ContigViz_Hexamer_continuous=ContigViz_Hexamer_continuous,ContigViz_Hexamer_digitized=ContigViz_Hexamer_digitized,ContigViz_FastK_continuous=ContigViz_FastK_continuous,ContigViz_FastK_digitized=ContigViz_FastK_digitized,ContigViz_Unique_15mers_continuous=ContigViz_Unique_15mers_continuous,ContigViz_Unique_15mers_digitized=ContigViz_Unique_15mers_digitized,ContigViz_Coverage_continuous=ContigViz_Coverage_continuous,ContigViz_Coverage_digitized=ContigViz_Coverage_digitized" --text-header {}'.format(cobiontid_contigviz_output_path, btk_dataset_folder)
        gpf.run_system_command(add_cobiontid_contigviz_results_command)



def main(assembly_fasta_path, dataset_folder, pipeline_run_folder, assembly_title, taxon_name, taxid, blastn_hits_path, uniprot_diamond_hits_path, nr_diamond_hits_path, mapped_reads_path, taxdump_path, threads, assembly_alias, dry_run_flag):

    #out_folder = pipeline_run_folder + "/collected_tables"

    if assembly_alias == "":
        assembly_alias = assembly_title

    if dry_run_flag == False:
        Path(dataset_folder).mkdir(parents=True, exist_ok=True)

    edited_assembly_title = assembly_title.replace(".", "_")
    edited_assembly_title = edited_assembly_title.replace(" ", "_")

    assembly_yaml_path = dataset_folder + "/" + edited_assembly_title + ".yaml"
    if dry_run_flag == False:
        create_assembly_yaml(assembly_yaml_path, assembly_alias, taxon_name)

    blobtools_create_command = "blobtools create --fasta {} --meta {} --taxid {} --taxdump {} {}".format(assembly_fasta_path, assembly_yaml_path, taxid, taxdump_path, dataset_folder)
    gpf.run_system_command(blobtools_create_command, dry_run=dry_run_flag)


    hits_file_paths = [blastn_hits_path, uniprot_diamond_hits_path, nr_diamond_hits_path]
    hits_file_paths = [n for n in hits_file_paths if os.path.isfile(n) is True and os.stat(n).st_size > 0]

    if len(hits_file_paths) > 0:
        add_hits_command = "blobtools add"
        for hits_file_path in hits_file_paths:
            add_hits_command += " --hits {}".format(hits_file_path)
        add_hits_command += " --taxrule bestsum --taxdump {} {}".format(taxdump_path, dataset_folder)
        gpf.run_system_command(add_hits_command, dry_run=dry_run_flag)


    if os.path.isfile(mapped_reads_path) is True and os.stat(mapped_reads_path).st_size > 0:
        add_cov_command = "blobtools add --cov {} --threads {} {}".format(mapped_reads_path, threads, dataset_folder)
        gpf.run_system_command(add_cov_command, dry_run=dry_run_flag)

    #export_table_command = "blobtools filter --table {}/btk_summary_table_basic.tsv {}".format(out_folder, dataset_folder)
    add_custom_variables_to_btk_dataset(pipeline_run_folder, dataset_folder)
    export_table_command = "blobtools filter --table {}/collected_tables/btk_summary_table_full.tsv {}".format(pipeline_run_folder, dataset_folder)

    gpf.run_system_command(export_table_command, dry_run=dry_run_flag)

    # json_file_path = run_folder + "/" + sample_id + ".json"
    #cobiontid_scaffs_json_to_csv(json_file_path, out_folder + "/cobiontid_markerscan.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("assembly_fasta_path", type=str, help="assembly_fasta_path")
    parser.add_argument("dataset_folder", type=str, help="Path for dataset folder")
    parser.add_argument("pipeline_run_folder", type=str, help="Folder where this pipeline is run pipeline")
    parser.add_argument("assembly_title", type=str, help="Assembly title")
    parser.add_argument("taxon_name", type=str, help="Taxon name")
    parser.add_argument("taxid", type=int, help="taxid")
    parser.add_argument("blastn_hits_path", type=str, help="Path to blastn hits file")
    parser.add_argument("uniprot_diamond_hits_path", type=str, help="Path to UNIPROT Diamond BLASTX hits file")
    parser.add_argument("nr_diamond_hits_path", type=str, help="Path to nr Diamond BLASTX hits file")
    parser.add_argument("mapped_reads_path", type=str, help="Path to the BAM file with mapped reads for coverage estimation")
    parser.add_argument("taxdump_path", type=str, help="Path to the directory with NCBI taxdump files")
    parser.add_argument("--threads", type=int, default=1, help="Number of CPU threads (default: 1)")
    parser.add_argument("--assembly_alias", type=str, default="", help="Assembly alias")
    parser.add_argument("--dry_run", dest="dry_run", action="store_true", help="Dry run (print commands without executing)")
    args = parser.parse_args()
    main(args.assembly_fasta_path, args.dataset_folder, args.pipeline_run_folder, args.assembly_title, args.taxon_name, args.taxid, args.blastn_hits_path, args.uniprot_diamond_hits_path, args.nr_diamond_hits_path, args.mapped_reads_path, args.taxdump_path, args.threads, args.assembly_alias,args.dry_run)
