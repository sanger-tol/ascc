#!/usr/bin/env python3
"""
Script for finding the average coverage of each scaffold in a SAMtools depth output file
"""

import sys
import general_purpose_functions as gpf

in_path = sys.argv[1]
in_data = gpf.ll(in_path)
scaffs_dict = dict()

for line in in_data:
    split_line = line.split()
    scaff_name = split_line[0]
    coverage = int(split_line[2])
    if scaff_name in scaffs_dict:
        scaffs_dict[scaff_name]["coverage_sum"] += coverage
        scaffs_dict[scaff_name]["scaff_len"] += 1
    else:
        scaffs_dict[scaff_name] = dict()
        scaffs_dict[scaff_name]["coverage_sum"] = coverage
        scaffs_dict[scaff_name]["scaff_len"] = 1

for scaff_name in scaffs_dict:
    average_coverage = scaffs_dict[scaff_name]["coverage_sum"] / scaffs_dict[scaff_name]["scaff_len"]
    print(scaff_name + "," + str(average_coverage))
