#!/usr/bin/env python
# This script was written by James Torrance (@jt8)
# Modified by Damon-Lee Pointon (@dp24 / @DLBPointon)

import csv
import re
import argparse

def main():

	parser = argparse.ArgumentParser(description='Generate text for ENA submission')
	parser.add_argument("--input", type=str, help='Input BED', default=None)
	parser.add_argument("--output", type=str, help='Output recommendation', default=None)
	parser.add_argument("-v", action='version', version='1.0')
	args = parser.parse_args()

	recommendation_handle = open(args.output, 'w')
	with open(args.input) as bed_handle:
		bed_csv_reader = csv.reader(bed_handle, delimiter='\t')
		for field_set in bed_csv_reader:
			scaffold = field_set[0]
			(length, percentage) = re.split(',', field_set[3])
			percentage = float(percentage)
			if percentage > 90:
				recommendation_handle.write(f'REMOVE\t{scaffold}\n')
			elif percentage > 50:
				recommendation_handle.write(f'Investigate\t{scaffold}\n')
	recommendation_handle.close()

if __name__ == "__main__":
		main()
