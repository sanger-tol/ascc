#!/usr/bin/env python3
"""
Script for finding the average coverage of each scaffold in a SAMtools depth output file.
originally written by Eerik Aunin (ea10)
re-wriiten by Yumi Sims (yy5)
"""

import sys
import argparse
import general_purpose_functions as gpf

# Define the version number
__version__ = '1.0.0'

def display_version():
    """Display the script's version."""
    print(f"{__version__}")

def process_data(in_data):
    """Process the input data and calculate scaffold coverage."""
    scaffs_dict = {}
    
    for line in in_data:
        split_line = line.split()
        scaff_name, coverage = split_line[0], int(split_line[2])
        
        scaffs_dict[scaff_name] = {
            "coverage_sum": scaffs_dict.get(scaff_name, {}).get("coverage_sum", 0) + coverage,
            "scaff_len": scaffs_dict.get(scaff_name, {}).get("scaff_len", 0) + 1
        }
    return scaffs_dict

def calculate_average_coverage(scaffs_dict):
    """Calculate average coverage for each scaffold."""
    results = []
    for scaff_name, data in scaffs_dict.items():
        average_coverage = data["coverage_sum"] / data["scaff_len"]
        results.append((scaff_name, average_coverage))
    return results

def main():
    """Main function to run the script."""
    parser = argparse.ArgumentParser(description="Calculate average coverage for each scaffold.")
    parser.add_argument('input_file', nargs='?', help='Path to the input SAMtools depth file')
    parser.add_argument('--version', action='store_true', help='Display the version')

    args = parser.parse_args()

    if args.version:
        display_version()
        return

    if not args.input_file:
        print("Error: Missing input file. Please provide an input file.")
        parser.print_help()
        sys.exit(1)

    in_data = gpf.ll(args.input_file)
    scaffs_dict = process_data(in_data)
    results = calculate_average_coverage(scaffs_dict)

    for scaff_name, average_coverage in results:
        print(f"{scaff_name}, {average_coverage:.2f}")

if __name__ == "__main__":
    main()