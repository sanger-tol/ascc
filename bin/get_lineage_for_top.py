#!/usr/bin/env python3
"""
Script for getting the lineage for BLAST hits, by combining a BLAST hits file with lineage information from NCBI's rankedlineage.dmp file
Actual format of the input row: qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore
"""

import csv
import argparse
import sys
import re


def parse_taxdump(taxdump_file):
    """
    Function to parse the rankedlineage.dmp file
    """
    taxid_to_lineage = {}
    with open(taxdump_file, "r") as f:
        for line in f:
            split_line = line.split("|")
            split_line = [n.strip() for n in split_line]
            taxid = split_line[0]
            lineage = split_line[1:]
            taxid_to_lineage[taxid] = lineage
    return taxid_to_lineage


def extract_taxid_from_accession(accession):
    """
    Function to extract taxid from accession if possible
    Returns a placeholder taxid if extraction fails
    """
    # This is a placeholder. In a real implementation, you would
    # query the NCBI database or use a mapping file to get the taxid.
    # For now, we'll use a placeholder taxid
    return "0"


def process_blast_results(blast_tsv, taxdump_data, output_csv, column_name_prefix):
    """
    Function to parse the BLAST TSV results and combine with taxdump data
    """
    with open(blast_tsv, "r") as blast_file, open(output_csv, "w", newline="") as output_file:
        blast_reader = csv.reader(blast_file, delimiter="\t")
        fieldnames = [
            "scaff", f"{column_name_prefix}blast_accession", f"{column_name_prefix}blast_score",
            f"{column_name_prefix}blast_taxid", f"{column_name_prefix}blast_species",
            f"{column_name_prefix}blast_genus", f"{column_name_prefix}blast_family",
            f"{column_name_prefix}blast_order", f"{column_name_prefix}blast_class",
            f"{column_name_prefix}blast_phylum", f"{column_name_prefix}blast_kingdom",
            f"{column_name_prefix}blast_superkingdom", f"{column_name_prefix}blast_no_rank",
            f"{column_name_prefix}blast_empty"
        ]
        output_writer = csv.DictWriter(output_file, fieldnames=fieldnames)
        output_writer.writeheader()

        for row in blast_reader:
            try:
                # Actual format: qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore
                if len(row) == 12:
                    qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = row
                    # Since taxid is missing, we'll use a placeholder
                    taxid = "0"
                else:
                    sys.stderr.write(f"Unexpected number of columns in row: {len(row)}\n")
                    continue

                # Create a default empty lineage with 10 elements (9 levels + empty)
                default_lineage = [""] * 10
                
                # If taxid is found in taxdump_data, use that lineage
                if taxid in taxdump_data:
                    lineage_data = taxdump_data[taxid]
                    # Copy available lineage data to our default lineage
                    for i in range(min(len(lineage_data), 9)):
                        default_lineage[i] = lineage_data[i]

                output_writer.writerow({
                    "scaff": qseqid,
                    f"{column_name_prefix}blast_accession": sseqid,
                    f"{column_name_prefix}blast_score": bitscore,
                    f"{column_name_prefix}blast_taxid": taxid,
                    f"{column_name_prefix}blast_species": default_lineage[0],
                    f"{column_name_prefix}blast_genus": default_lineage[1],
                    f"{column_name_prefix}blast_family": default_lineage[2],
                    f"{column_name_prefix}blast_order": default_lineage[3],
                    f"{column_name_prefix}blast_class": default_lineage[4],
                    f"{column_name_prefix}blast_phylum": default_lineage[5],
                    f"{column_name_prefix}blast_kingdom": default_lineage[6],
                    f"{column_name_prefix}blast_superkingdom": default_lineage[7],
                    f"{column_name_prefix}blast_no_rank": default_lineage[8],
                    f"{column_name_prefix}blast_empty": default_lineage[9],
                })
            except Exception as e:
                sys.stderr.write(f"Error processing row: {e}\n")
                continue

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process BLAST results and map taxonomic lineages.")
    parser.add_argument("--blast_tsv", required=True, help="Input BLAST TSV file")
    parser.add_argument("--taxdump", required=True, help="Input rankedlineage.dmp file")
    parser.add_argument("--output_csv", required=True, help="Output CSV file")
    parser.add_argument("--column_name_prefix", default="nt_", help="Column name prefix (default: nt_)")

    args = parser.parse_args()

    taxdump_data = parse_taxdump(args.taxdump)
    process_blast_results(args.blast_tsv, taxdump_data, args.output_csv, args.column_name_prefix)

    sys.stderr.write(f"Output written to {args.output_csv}\n")
