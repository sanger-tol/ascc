import sys
import re
import os


class BedTools:
    # bedtools = '/software/grit/tools/bedtools-2.29.0/bin/bedtools'
    bedtools = "bedtools"

    def sort_and_merge_bed_file(self, bed_file):
        sorted_bed_file = re.sub("\.bed$", ".sorted.bed", bed_file)
        merged_bed_file = re.sub("\.bed$", ".merged.bed", bed_file)
        sort_command = self.bedtools + " sort -i " + bed_file + " > " + sorted_bed_file
        print(sort_command)
        os.system(sort_command)
        merge_command = self.bedtools + " merge -i " + sorted_bed_file + " > " + merged_bed_file
        print(merge_command)
        os.system(merge_command)
        return merged_bed_file

    def merge_bed_b_into_a(self, bed_file_1, bed_file_2):
        # We're merging file_2 into file_1
        os.system("cat " + bed_file_2 + " >> " + bed_file_1)
        merged_bed_file = sort_and_merge_bed_file(bed_file_1)
        os.system("mv " + merged_bed_file + " " + bed_file_1)

    def subtract_b_from_a(self, bed_file_1, bed_file_2):
        # We subtract file 2 from file 1
        subtracted_bed_file = re.sub("\.bed$", ".subtracted.bed", bed_file_1)
        subtract_command = (
            self.bedtools + " subtract -a " + bed_file_1 + " -b " + bed_file_2 + " > " + subtracted_bed_file
        )
        print(subtract_command)
        os.system(subtract_command)
        return subtracted_bed_file

    def coverage_for_bed_file(self, bed_file):
        coverage = 0

        with open(bed_file, "r") as bed_handle:
            for line in bed_handle:
                line = line.rstrip()
                fields = re.split("\s+", line)
                coverage += int(fields[2]) - int(fields[1])

        return coverage

    def coverage_for_bed_file_by_scaffold(self, bed_file):
        coverage_for_scaffold = {}

        with open(bed_file, "r") as bed_handle:
            for line in bed_handle:
                line = line.rstrip()
                fields = re.split("\s+", line)
                if fields[0] not in coverage_for_scaffold:
                    coverage_for_scaffold[fields[0]] = 0
                coverage_for_scaffold[fields[0]] += int(fields[2]) - int(fields[1])

        return coverage_for_scaffold

    def coords_to_bed(self, coord_list_for_sequence, bed_file):
        bed_handle = open(bed_file, "w")
        for sequence in coord_list_for_sequence:
            for coord_pair in coord_list_for_sequence[sequence]:
                bed_handle.write("\t".join([sequence, str(coord_pair[0] - 1), str(coord_pair[1])]) + "\n")
        bed_handle.close()

    def bed_to_coords(self, bed_file):
        coord_list_for_sequence = {}

        with open(bed_file, "r") as bed_handle:
            for line in bed_handle:
                line = line.rstrip()
                fields = re.split("\s+", line)
                if fields[0] not in coord_list_for_sequence:
                    coord_list_for_sequence[fields[0]] = []
                coord_list_for_sequence[fields[0]].append([int(fields[1]) + 1, int(fields[2])])

        return coord_list_for_sequence
