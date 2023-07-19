#!/usr/bin/env python3

"""
Script for retrieving lineage information for BLAST hits

Written by Eerik Aunin @eeaunin

Adapted by Damon-Lee Pointon @DLBPointon
"""

import sys
import general_purpose_functions as gpf
import argparse
import os

def extract_accessions_from_blast_output(blast_output_path, accessions_set_out_path):
    """
    Input: 1) path to BLAST output file (format 6), 2) path for output file
    Output: a set of BLAST hit accession numbers from the input file, without accession version number.
    This is written to a text file, each element on a separate row.
    """
    in_data = gpf.ll(blast_output_path)

    accessions_list = list()
    for line in in_data:
        split_line = line.split()
        accession = split_line[1]
        split_accession = accession.split(".")
        truncated_accession = split_accession[0]
        accessions_list.append(truncated_accession)

    accessions_set = set(accessions_list)
    gpf.export_list_as_line_break_separated_file(list(accessions_set), accessions_set_out_path)


def get_taxid_by_accession(accessions_set_out_path, accession2taxid_folder):
    """
    Input: 1) path to output file of extract_accessions_from_blast_output function, 2) path to accession2taxid files downloaded from NCBI FTP
    Output: dictionary where keys are accession numbers (without version number) and values are corresponding taxIDs
    """
    blast_accessions = gpf.l(accessions_set_out_path)
    accession2taxid_files = os.listdir(accession2taxid_folder)
    accession2taxid_files = [n for n in accession2taxid_files if n.endswith(".accession2taxid")]
    if len(accession2taxid_files) == 0:
        sys.stderr.write(f"No *.accession2taxid files were found in the expected location ({accession2taxid_folder})\n")
        sys.exit(1)
    out_dict = dict()
    for accession2taxid_file in accession2taxid_files:
        accession2taxid_data = gpf.ll(f"{accession2taxid_folder}/{accession2taxid_file}")
        for line in accession2taxid_data:
            split_line = line.split("\t")
            assert len(split_line) == 4
            accession = split_line[0]
            taxid = split_line[2]

            if accession in blast_accessions:
                out_dict[accession] = taxid
    return out_dict


def get_accession_to_lineage_dict(accession2taxid_dict, lineages_path):
    """
    Input:  1) a dictionary where keys are accession numbers (without version number) and values are corresponding taxIDs
            2) path to a lineages file (rankedlineage.dmp)
    Output: a dictionary where keys are taxIDs and the values are the corresponding lineage data
    """
    lineages_data = gpf.ll(lineages_path)

    taxid2lineage_dict = dict()
    taxid_values = list(accession2taxid_dict.values())
    taxid_values = list(set(taxid_values))
    for line in lineages_data:
        split_line = line.split("|")
        split_line = [n.strip() for n in split_line]
        taxid = split_line[0]
        if taxid in taxid_values:
            lineage = ",".join(split_line[1:len(split_line)])
            taxid2lineage_dict[taxid] = lineage
    return taxid2lineage_dict


def parse_blast_results(blast_output_path):
    """
    Input: path to BLAST output file (format 6)
    Output: a dictionary where the keys are BLAST query contig names. The values are dictionaries where the keys are BLAST hit accession numbers and the
        values are BLAST scores of the corresponding hits
    """
    results_dict = dict()
    with open(blast_output_path, "r") as in_file:
        for line in in_file:
            old_top_score = None
            line = line.rstrip()
            split_line = line.split()
            query_contig = split_line[0]
            blast_hit_accession = split_line[1]
            blast_hit_score = float(split_line[11])
            if query_contig not in results_dict:
                results_dict[query_contig] = {blast_hit_accession: blast_hit_score}
            else:
                if blast_hit_accession in results_dict[query_contig]:
                    old_top_score = results_dict[query_contig][blast_hit_accession]

                    if blast_hit_score > old_top_score:
                        results_dict[query_contig][blast_hit_accession] = blast_hit_score
                else:
                    results_dict[query_contig][blast_hit_accession] = blast_hit_score
    return results_dict


def main(blast_output_path, out_folder, accession2taxid_folder, lineages_path, column_name_prefix):
    accessions_set_out_path = out_folder + "/BLAST_hit_accessions.txt"
    extract_accessions_from_blast_output(blast_output_path, accessions_set_out_path)
    accession2taxid_dict = get_taxid_by_accession(accessions_set_out_path, accession2taxid_folder)
    taxid2lineage_dict = get_accession_to_lineage_dict(accession2taxid_dict, lineages_path)
    blast_results_dict = parse_blast_results(blast_output_path)

    out_path = out_folder + "/" + "BLAST_results_with_lineage.csv"
    out_list = list()
    if column_name_prefix != "":
        column_name_prefix += "_"
    out_list.append("scaff,{0}blast_accession,{0}blast_score,{0}blast_taxid,{0}blast_species,{0}blast_genus,{0}blast_family,{0}blast_order,{0}blast_class,{0}blast_phylum,{0}blast_kingdom,{0}blast_superkingdom,{0}blast_no_rank,{0}blast_empty".format(column_name_prefix))
    for key in blast_results_dict:
        contig_results_dict = blast_results_dict[key]
        for result_key in contig_results_dict:
            accession = str(result_key)
            score = str(contig_results_dict[result_key])
            lineage = "NA"
            taxid = "NA"
            truncated_accesion = accession
            if "." in accession:
                truncated_accesion = accession.split(".")[0]
            if truncated_accesion in accession2taxid_dict:
                taxid = accession2taxid_dict[truncated_accesion]
                if taxid in taxid2lineage_dict:
                    lineage = taxid2lineage_dict[taxid]
            out_line = ",".join([str(key), accession, score, taxid, lineage])
            out_list.append(out_line)

    gpf.export_list_as_line_break_separated_file(out_list, out_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-v", action='version', version='1.0')
    parser.add_argument("blast_output_path", type=str, help="Path to BLAST output file (format 6)")
    parser.add_argument("out_folder", type=str, help="Folder for output files")
    parser.add_argument("accession2taxid_folder", type=str, help="Path to folder with *.accession2taxid files, downloaded from NCBI FTP")
    parser.add_argument("lineages_path", type=str, help="Path NCBI rankedlineage.dmp file, downloaded from NCBI FTP")
    parser.add_argument("--column_name_prefix", type=str, help="Prefix that will be added to column names in the output CSV file. Default: empty string", default="")
    args = parser.parse_args()
    main(args.blast_output_path, args.out_folder, args.accession2taxid_folder, args.lineages_path, args.column_name_prefix)
