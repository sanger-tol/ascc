#!/usr/bin/env python3
"""
This script is based on an AWK script from NCBI: https://ftp.ncbi.nlm.nih.gov/pub/kitts/VSlistTo1HitPerLine.awk, orginally by Paul Kitts.
Converted to Python by Eerik Aunin.

This script converts the VecScreen text list output to one line giving the coordinates for each vector segment in the format:
VecScreen_Category   ID_string   start_position   end_position
The default is to report Strong, Moderate, and Weak matches and also segments of Suspect Origin. Reporting of any category can be suppressed by including 
--skip_reporting_suspect_hits, --skip_reporting_weak_hits, --skip_reporting_moderate_hits or --skip_reporting_strong_hits on the command line. 
"No hits" will be reported for any Query sequence that had no matches in any of the selected categories, unless --skip_reporting_no_hits is included on the command line.
VecScreen errors will be reported unless --skip_reporting_errors is included on the command line.
Usage:
VSlistTo1HitPerLine.py [options] vecscreen_output_file
"""

import re
import argparse

def main(args):
    hits_to_report = ""
    ID = ""
    hits = 0
    error_found = 0

    with open(args.vecscreen_output_file, "r") as file:
        for line in file:
            line = line.strip()
            Fld = line.split(" ")

            if not re.match(r'^[0-9 \t]+$', line):
                hits_to_report = ""

            if hits_to_report:
                print(f"VecScreen_{hits_to_report.ljust(8)}\t{ID}\t{line}")
                hits += 1
                continue
            
            if line.startswith(">Vector "):
                if ID != "":
                    if error_found:
                        print(f"VecScreen_ERROR   \t{ID}")
                    elif hits == 0 and args.skip_reporting_no_hits is False:
                        print(f"VecScreen_No_Hits \t{ID}")
                ID = Fld[1]
                hits = 0
                error_found = 0
                continue

            if args.skip_reporting_strong_hits is False and line.startswith("Strong"):
                hits_to_report = "Strong"
                continue

            if args.skip_reporting_moderate_hits is False and line.startswith("Moderate"):
                hits_to_report = "Moderate"
                continue

            if args.skip_reporting_weak_hits is False and line.startswith("Weak"):
                hits_to_report = "Weak"
                continue

            if args.skip_reporting_suspect_hits is False and line.startswith("Suspect"):
                hits_to_report = "Suspect"
                continue

            if args.skip_reporting_errors is False:
                if line.startswith("ERROR") or line.startswith("WARNING"):
                    error_found += 1

    if ID != "":
        if error_found > 0:
            print(f"VecScreen_ERROR   \t{ID}")
        elif hits == 0 and args.skip_reporting_no_hits is False:
            print(f"VecScreen_No_Hits \t{ID}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reformatting VecScreen's output")
    parser.add_argument("vecscreen_output_file", type=str, help="Path to a raw output file from NCBI VecScreen (from a run with the -f3 flag)", default=None)
    parser.add_argument("--skip_reporting_strong_hits", action="store_true", help="Skip reporting strong hits")
    parser.add_argument("--skip_reporting_moderate_hits", action="store_true", help="Skip reporting moderate hits")
    parser.add_argument("--skip_reporting_weak_hits", action="store_true", help="Skip reporting weak hits")
    parser.add_argument("--skip_reporting_suspect_hits", action="store_true", help="Skip reporting hits of suspect origin")
    parser.add_argument("--skip_reporting_no_hits", action="store_true", help="Skip reporting no-hits")
    parser.add_argument("--skip_reporting_errors", action="store_true", help="Skip reporting errors")
    parser.add_argument("-v", "--version", action="version", version="1.0")
    args = parser.parse_args()
    main(args)