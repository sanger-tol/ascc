#!/usr/bin/env python3
"""
Script for finding the average coverage of each scaffold in a SAMtools depth output file.
originally written by Eerik Aunin (ea10)
refactored by Yumi Sims (yy5)
"""

import sys
import argparse
import general_purpose_functions as gpf


def process_data(in_data):
    """Process the input data and calculate scaffold coverage."""
    scaffs_dict = {}

    for line in in_data:
        split_line = line.split()
        scaff_name, coverage = split_line[0], int(split_line[2])

        scaffs_dict[scaff_name] = {
            "coverage_sum": scaffs_dict.get(scaff_name, {}).get("coverage_sum", 0)
            + coverage,
            "scaff_len": scaffs_dict.get(scaff_name, {}).get("scaff_len", 0) + 1,
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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("in_path", type=str, help="Path to SAMtools depth output file")
    parser.add_argument("-v", "--version", action="version", version="1.0")
    args = parser.parse_args()

    try:
        in_data = gpf.ll(args.in_path)
        scaffs_dict = process_data(in_data)
        results = calculate_average_coverage(scaffs_dict)

        for scaff_name, average_coverage in results:
            print(f"{scaff_name}, {average_coverage:.2f}")
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
