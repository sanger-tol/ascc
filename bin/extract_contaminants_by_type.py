#!/usr/bin/env python
"""
Script from James Torrance (jt8)
Minor modifications by Eerik Aunin (ea10) and Damon-Lee Pointon (@dp24/@DLBPointon)
"""

import os
import re
import BedTools
from Bio import SeqIO
import gzip
import argparse
import general_purpose_functions as gpf


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("contamination_files", nargs="+")
    parser.add_argument("--assembly_file")
    parser.add_argument("--extracts_dir", default="./")
    parser.add_argument("-v", action="version", version="1.0")
    args = parser.parse_args()

    for contamination_file in args.contamination_files:
        print("Getting data for " + contamination_file)

        if args.assembly_file == None:
            assembly_file = re.sub("\.contamination", ".fa.gz", contamination_file)
        else:
            assembly_file = args.assembly_file
        assembly_name = "assembly"

        # If at first you don't succeed, try again with a .fasta suffix:
        if not os.path.isfile(assembly_file):
            assembly_file = re.sub("\.contamination", ".fasta.gz", contamination_file)

        if not os.path.isfile(assembly_file):
            print(assembly_file + " does not exist")
        else:
            coord_list_for_section_and_sequence = {
                "ALL": {},
            }

            section_name = "NO_SECTION"

            contamination_handle = open(contamination_file, "r")
            for line in contamination_handle:
                section_name_match = re.search("^\=+\s*(.+)\s+\=", line)
                if section_name_match:
                    section_name = section_name_match.group(1)
                    if section_name in coord_list_for_section_and_sequence:
                        exit("Duplicate section name: " + section_name + "\n")

                if not (re.search("^[#=]", line)):
                    fields = line.split("\t")
                    if len(fields) > 7:
                        if section_name not in coord_list_for_section_and_sequence:
                            coord_list_for_section_and_sequence[section_name] = {}
                        if fields[0] not in coord_list_for_section_and_sequence[section_name]:
                            coord_list_for_section_and_sequence[section_name][fields[0]] = []
                        if fields[0] not in coord_list_for_section_and_sequence["ALL"]:
                            coord_list_for_section_and_sequence["ALL"][fields[0]] = []

                        coords = [int(fields[6]), int(fields[7])]
                        coords = sorted(coords)
                        coord_list_for_section_and_sequence[section_name][fields[0]].append(coords)
                        coord_list_for_section_and_sequence["ALL"][fields[0]].append(coords)

            # margins = [0,10000]
            margins = [0]

            # Write BED file
            for section_name in coord_list_for_section_and_sequence:
                # Only print sections that have data- but always produce an "ALL" file
                if len(coord_list_for_section_and_sequence[section_name]) > 0 or section_name == "ALL":
                    output_section_name = section_name
                    output_section_name = re.sub("(\s|\:)", "_", section_name)
                    bed_file = args.extracts_dir + assembly_name + "." + output_section_name + ".bed"
                    bedtools = BedTools.BedTools()
                    bedtools.coords_to_bed(coord_list_for_section_and_sequence[section_name], bed_file)
                    merged_bed_file = bedtools.sort_and_merge_bed_file(bed_file)
                    merged_coord_list_for_sequence = bedtools.bed_to_coords(merged_bed_file)

                    # Get lengths
                    length_file = args.extracts_dir + assembly_name + ".lengths"

                    # Replaced the exonerate fastalength binary with a Python script (EA)
                    fastalength(assembly_file, length_file)
                    length_for_sequence = parse_fastalength_file(length_file)

                    coverage_file_base_name = args.extracts_dir + assembly_name + "." + output_section_name
                    write_coverage_file(coverage_file_base_name, merged_bed_file, length_for_sequence, bedtools)

                    for margin in margins:
                        output_file = args.extracts_dir + assembly_name + "." + output_section_name + ".extracts.2.fa"
                        output_handle = open(output_file, "w")

                        handle = None

                        if re.search("gz$", assembly_file):
                            handle = gzip.open(assembly_file, "rt")
                        else:
                            handle = open(assembly_file, "rt")

                        for record in SeqIO.parse(handle, "fasta"):
                            if record.id in merged_coord_list_for_sequence:
                                for coord_pair in merged_coord_list_for_sequence[record.id]:
                                    extracted_sequence = extract_sequence(record, coord_pair[0], coord_pair[1], margin)
                                    SeqIO.write([extracted_sequence], output_handle, "fasta")
                        handle.close()


def write_coverage_file(coverage_file_base_name, merged_bed_file, length_for_sequence, bedtools):
    # Record coverage
    merged_coord_list_for_sequence = bedtools.bed_to_coords(merged_bed_file)
    coverage_for_sequence = bedtools.coverage_for_bed_file_by_scaffold(merged_bed_file)

    percentage_coverage_for_sequence = {}
    for sequence in coverage_for_sequence:
        if sequence in length_for_sequence:
            percentage_coverage_for_sequence[sequence] = (
                coverage_for_sequence[sequence] / length_for_sequence[sequence] * 100
            )
        else:
            # Handle missing sequence length information
            percentage_coverage_for_sequence[sequence] = 0

    coverage_threshold = 1

    # Take the stem name for the file, and print to various files, with c coverage threshold, without, and per line

    filtered_coverage_scaffold_file = coverage_file_base_name + ".filtered_scaffold_coverage.bed"
    unfiltered_coverage_scaffold_file = coverage_file_base_name + ".unfiltered_scaffold_coverage.bed"

    filtered_coverage_scaffold_handle = open(filtered_coverage_scaffold_file, "w")
    unfiltered_coverage_scaffold_handle = open(unfiltered_coverage_scaffold_file, "w")

    bed_with_coverage_format = "{0}\t{1}\t{2}\t{3:d},{4:.3f}\n"

    for sequence in sorted(
        percentage_coverage_for_sequence.keys(), key=lambda x: percentage_coverage_for_sequence[x], reverse=True
    ):
        if sequence not in length_for_sequence:
            # NB Query names including .suffixes may have their suffixes lost in the BLAST stage
            # We're not usually including such names, so this isn't an issue at the moment
            exit("Cannot find length for " + sequence)
        unfiltered_coverage_scaffold_handle.write(
            bed_with_coverage_format.format(
                sequence,
                0,
                length_for_sequence[sequence],
                coverage_for_sequence[sequence],
                percentage_coverage_for_sequence[sequence],
            )
        )
        if percentage_coverage_for_sequence[sequence] >= coverage_threshold:
            filtered_coverage_scaffold_handle.write(
                bed_with_coverage_format.format(
                    sequence,
                    0,
                    length_for_sequence[sequence],
                    coverage_for_sequence[sequence],
                    percentage_coverage_for_sequence[sequence],
                )
            )
    unfiltered_coverage_scaffold_handle.close()
    filtered_coverage_scaffold_handle.close()

    filtered_coverage_region_file = coverage_file_base_name + ".filtered_region_coverage.bed"
    unfiltered_coverage_region_file = coverage_file_base_name + ".unfiltered_region_coverage.bed"

    filtered_coverage_region_handle = open(filtered_coverage_region_file, "w")
    unfiltered_coverage_region_handle = open(unfiltered_coverage_region_file, "w")

    for sequence in merged_coord_list_for_sequence:
        for coord_pair in merged_coord_list_for_sequence[sequence]:
            unfiltered_coverage_region_handle.write(
                bed_with_coverage_format.format(
                    sequence,
                    str(coord_pair[0] - 1),
                    str(coord_pair[1]),
                    coverage_for_sequence[sequence],
                    percentage_coverage_for_sequence[sequence],
                )
            )
            if percentage_coverage_for_sequence[sequence] >= coverage_threshold:
                filtered_coverage_region_handle.write(
                    bed_with_coverage_format.format(
                        sequence,
                        str(coord_pair[0] - 1),
                        str(coord_pair[1]),
                        coverage_for_sequence[sequence],
                        percentage_coverage_for_sequence[sequence],
                    )
                )


def fastalength(input_path, length_file):
    fasta_data = gpf.read_fasta_in_chunks(input_path)
    with open(length_file, "w") as length_file_handle:
        for header, seq in fasta_data:
            if header is not None:  # Check if header is not None
                if " " in header:
                    header = header.replace(" ", "_")
                result = f"{len(seq)} {header}"
                print(result)
                length_file_handle.write(result + "\n")


# Parse fastalength file
def parse_fastalength_file(length_file):
    length_for_sequence = {}

    with open(length_file, "r") as length_handle:
        for line in length_handle:
            line = line.rstrip()
            fields = re.split("\s+", line)
            length_for_sequence[fields[1]] = int(fields[0])

    return length_for_sequence


# Extract sequence from assembly
def extract_sequence(seq_record, start, end, margin):
    # Extend the seq-record by the margin
    # First add to end
    if end + margin < len(seq_record.seq):
        end += margin
    else:
        end = len(seq_record.seq)
    # Deduct from start
    if start - margin > 1:
        start -= margin
    else:
        start = 1

    extracted_slice = seq_record[start:end]
    extracted_slice.id = seq_record.id + ":" + str(start) + "-" + str(end)
    extracted_slice.name = extracted_slice.id

    return extracted_slice


if __name__ == "__main__":
    main()
