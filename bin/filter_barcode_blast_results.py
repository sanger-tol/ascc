#!/usr/bin/env python3
"""
Script from James Torrance for filtering the results of running BLAST
with assembly FASTA file against PacBio barcode sequences.

Originally written James Torrance (@jt8)

Modified by Damon-Lee Pointon (@dp24 / @DLBPointon)
"""

import re
import gzip
from Bio import SeqIO
import argparse

def main():

        parser = argparse.ArgumentParser(description='Identify barcodes for removal')
        parser.add_argument("--input", help='FASTA input')
        parser.add_argument("--blast", help='BLAST data')
        parser.add_argument("--output", help='FASTA output')
        parser.add_argument("--barcode", help='Barcode to filter for')
        parser.add_argument("-v", action='version', version='1.0')
        args = parser.parse_args()

        barcode_location_file = args.output
        if args.output == None:
                barcode_location_file = args.blast + '.debarcoded'

        barcode_location_handle = open(barcode_location_file, 'w')

        n_regions_for_record = {}
        length_for_record = {}

        fasta_input_handle = None
        if re.search('\.gz$', args.input):
                fasta_input_handle = gzip.open(args.input, "rt")
        else:
                fasta_input_handle = open(args.input, "rt")

        for record in SeqIO.parse(fasta_input_handle, "fasta"):
                # barcode_location_handle.write(record.id + '\n')
                n_regions_for_record[record.id] = []
                length_for_record[record.id] = len(str(record.seq))
                n_iterator = re.finditer('[Nn]+', str(record.seq))
                for n_instance in n_iterator:
                        # barcode_location_handle.write(record.id + '\t' + str(n_instance.start(0)+1) + '\t' + str(n_instance.end(0)) + '\n')
                        n_regions_for_record[record.id].append([n_instance.start(0)+1, n_instance.end(0)])

        fasta_input_handle.close()
        with open(args.blast, "r") as blast_input_handle:
                for line in blast_input_handle:
                        if not re.search('^#', line):
                                #barcode_location_handle.write(line)
                                blast_fields = re.split('\s+', line)
                                blast_id = blast_fields[0]
                                blast_barcode_id = blast_fields[1]
                                blast_percentage = float(blast_fields[2])
                                blast_length = int(blast_fields[3])
                                blast_start = int(blast_fields[6])
                                blast_end = int(blast_fields[7])

                                # Does this have the expected barcode?
                                if args.barcode == None or re.match(args.barcode, blast_barcode_id):
                                        # Is this within N of a start or end?
                                        starts = [1]
                                        ends = [length_for_record[blast_id]]
                                        for n_region in n_regions_for_record[blast_id]:
                                                # barcode_location_handle.write('Adding N ' + str(n_region[0]) + '-' + str(n_region[1]) + '\n')
                                                starts.append(n_region[1]+1)
                                                ends.append(n_region[0]-1)
                                        threshold = 10
                                        terminal_flag = False
                                        for start in starts:
                                                if blast_start >= start and blast_start - start <= threshold:
                                                        terminal_flag = True
                                                        blast_start = start
                                        for end in ends:
                                                if blast_end <= end and end - blast_end <= threshold:
                                                        terminal_flag = True
                                                        blast_end = end
                                        if terminal_flag or (blast_length >= 17 and blast_percentage >= 100):
                                                barcode_location_handle.write('\t'.join(['CONTAMINANT',blast_id, str(blast_start), str(blast_end)]) + '\n')

        barcode_location_handle.close()
        # Generate N-coords
        # Parse

if __name__ == "__main__":
                main()

# MERGE COORDS
# COLLECT ALL CONTAMINATION/FASTA FILES FROM DIRECTORIES (RE-USE OLD CODE)
