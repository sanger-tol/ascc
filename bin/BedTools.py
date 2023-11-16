import sys
import re
import os
import pybedtools
from collections import defaultdict
from itertools import groupby


"""
Script from James Torrance (jt8) and Eerik Aunin (ea10)
refactored by Yumi sims (yy5)
"""


class BedTools:
    def __init__(self):
        pass

    def sort_and_merge_bed_file(self, input_bed_file):
        bed = pybedtools.BedTool(input_bed_file)
        # Get the output merged BED file name without the file extension
        output_merged_bed_file = os.path.splitext(os.path.basename(input_bed_file))[0] + ".merged.bed"
        # Sort the BED file
        sorted_bed = bed.sort()
        # Merge the sorted BED file
        merged_bed = sorted_bed.merge()
        # Save the merged BED file
        merged_bed.saveas(output_merged_bed_file)
        return output_merged_bed_file

    def merge_bed_b_into_a(self, bed_file_1, bed_file_2):
        # We're merging file_2 into file_1
        os.system("cat " + bed_file_2 + " >> " + bed_file_1)
        # Call the sort_and_merge_bed_file method to sort and merge the resulting file
        merged_bed_file = self.sort_and_merge_bed_file(bed_file_1)
        os.system("mv " + merged_bed_file + " " + bed_file_1)

    def subtract_b_from_a(self, bed_file_1, bed_file_2):
        if not os.path.exists(bed_file_1) or not os.path.exists(bed_file_2):
            raise FileNotFoundError("One or both of the input BED files do not exist.")
        output_file_name = os.path.splitext(os.path.basename(bed_file_1))[0] + ".subtracted.bed"

        # Load the BED files
        bed1 = pybedtools.BedTool(bed_file_1)
        bed2 = pybedtools.BedTool(bed_file_2)
        # Subtract bed_file_2 from bed_file_1
        subtracted_bed = bed1.subtract(bed2)
        subtracted_bed.saveas(output_file_name)
        return output_file_name

    def coverage_for_bed_file(self, bed_file):
        with open(bed_file, "r") as bed_handle:
            lines = [line.strip() for line in bed_handle]
            coverage = sum(int(fields[2]) - int(fields[1]) for line in lines for fields in [re.split("\s+", line)])
        return coverage

    def coverage_for_bed_file_by_scaffold(self, bed_file):
        with open(bed_file, "r") as bed_handle:
            fields_list = [re.split("\s+", line.strip()) for line in bed_handle]
        coverage_for_scaffold = {}
        for chrom, interval_list in groupby(fields_list, key=lambda x: x[0]):
            total_coverage = 0
            for _, start_str, end_str in interval_list:
                try:
                    start, end = int(start_str), int(end_str)
                    total_coverage += end - start
                except ValueError:
                    pass
            coverage_for_scaffold[chrom] = total_coverage
        return coverage_for_scaffold

    def coords_to_bed(self, coord_list_for_sequence, bed_file):
        with open(bed_file, "w") as bed_handle:
            lines = [
                f"{sequence}\t{start - 1}\t{end}\n"
                for sequence, coord_pairs in coord_list_for_sequence.items()
                for start, end in coord_pairs
            ]
            bed_handle.writelines(lines)

    def bed_to_coords(self, bed_file):
        coord_list_for_sequence = defaultdict(list)

        with open(bed_file, "r") as bed_handle:
            for line in bed_handle:
                fields = line.split()
                if len(fields) >= 3:
                    try:
                        seq_name, start, end = fields[0], int(fields[1]), int(fields[2])
                        coord_list_for_sequence[seq_name].append([start + 1, end])
                    except (ValueError, IndexError):
                        # Handle invalid lines or non-numeric fields
                        pass

        return dict(coord_list_for_sequence)
