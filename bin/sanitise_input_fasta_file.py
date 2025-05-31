#!/usr/bin/env python3
VERSION = "1.3.0"
DESCRIPTION = f"""
---
Script for sanitising FASTA headers and sequences:
- Shortens headers by splitting by whitespace and keeping only the first element
- Replaces problematic characters in headers (commas, spaces, etc.) with underscores
- Converts sequences to uppercase and replaces non-ATGC bases with N
Version: {VERSION}
---

Written by Eerik Aunin (ea10)
Modified by Damon-Lee Pointon (@dp24/@DLBPointon)
Further modified by Eerik Aunin (ea10)
Additional logging functionality added

"""

# MIT License
#
# Copyright (c) 2020-2022 Genome Research Ltd.
#
# Author: Eerik Aunin (eeaunin@gmail.com)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import general_purpose_functions as gpf
import argparse
import textwrap
import sys
import tempfile
import re
import json
import os


def is_all_n_sequence(seq):
    """Return True if sequence consists entirely of N's."""
    return all(base == "N" for base in seq.strip().upper())


def sanitise_sequence(seq):
    """Convert sequence to uppercase and replace any non-ATGC bases with N."""
    seq = seq.upper()
    return re.sub(r"[^ATGC]", "N", seq)


def sanitise_header(header):
    """Replace problematic characters in FASTA headers with underscores."""
    # Remove the '>' character if present at the start
    if header.startswith(">"):
        header = header[1:]

    # Replace problematic characters with underscores
    sanitised = re.sub(r"[,;\s|:]", "_", header)

    # Add back the '>' character
    return ">" + sanitised


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="sanitise_input_fasta_file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(DESCRIPTION),
    )
    parser.add_argument("fasta_path", type=str, help="Path to input FASTA file")
    parser.add_argument(
        "--delimiter",
        type=str,
        help="Delimiter string for splitting FASTA headers. Default: any whitespace character",
        default="",
    )
    parser.add_argument(
        "--allow_duplicate_headers", dest="allow_duplicate_headers", action="store_true"
    )
    parser.add_argument(
        "--keep_n_sequences",
        action="store_true",
        help="Keep sequences that are all Ns (default: False)",
    )
    parser.add_argument(
        "--log_file",
        type=str,
        help="Path to output log file in JSON format",
        default=None,
    )
    parser.add_argument(
        "--max_detailed_changes",
        type=int,
        help="Maximum number of detailed changes to report (0 for none)",
        default=100,
    )
    parser.add_argument(
        "--max_header_len",
        type=int,
        help="Maximum allowed length for FASTA headers (excluding '>'). If exceeded, script will exit with error.",
        default=0,  # Default 0 means no check
    )
    parser.add_argument("-v", "--version", action="version", version=VERSION)
    return parser.parse_args(argv)


def main(
    fasta_path,
    delimiter,
    allow_duplicate_headers,
    keep_n_sequences=False,
    log_file=None,
    max_detailed_changes=100,
    max_header_len=0,
):
    # Initialize log statistics
    log_stats = {
        "input_file": os.path.basename(fasta_path),
        "total_sequences": 0,
        "headers_shortened": 0,
        "headers_with_problematic_chars": 0,
        "sequences_with_non_atgc": 0,
        "non_atgc_bases_replaced": 0,
        "all_n_sequences_skipped": 0,
        "duplicate_headers_detected": 0,
        "headers_exceeding_max_length": 0,
        "has_issues": False,
        "detailed_changes": [],
    }

    with tempfile.TemporaryDirectory() as tmp_dir:
        input_file = fasta_path
        if fasta_path.endswith(".gz") or fasta_path.endswith('.gz"'):
            input_file = "{}/input_file.fa".format(tmp_dir)
            gpf.run_system_command("gunzip -c {} > {}".format(fasta_path, input_file))

        headers_list = list()
        headers_with_commas = 0
        in_data = gpf.ll(input_file)

        current_header = None
        current_sequence = []
        original_headers = {}  # Store original headers for logging

        def process_sequence():
            if current_header and current_sequence:
                sequence = "".join(current_sequence)
                log_stats["total_sequences"] += 1

                # Check if sequence is all N's
                if not keep_n_sequences and is_all_n_sequence(sequence):
                    log_stats["all_n_sequences_skipped"] += 1
                    log_stats["has_issues"] = True
                    log_stats["detailed_changes"].append(
                        {
                            "header": current_header[1:],
                            "change_type": "skipped_all_n_sequence",
                            "details": "Sequence consisted entirely of N's and was skipped",
                        }
                    )
                    sys.stderr.write(
                        "Skipping all-N sequence: {}\n".format(
                            current_header[1:].strip()
                        )
                    )
                else:
                    print(current_header)
                    print(sequence)

        for line in in_data:
            if line.startswith(">"):
                # Process previous sequence if it exists
                process_sequence()

                # Start new sequence
                original_header = line.strip()

                # Store original header for logging
                shortened_header = None
                if delimiter == "":
                    shortened_header = original_header.split()[0]
                else:
                    shortened_header = original_header.split(delimiter)[0]

                # Maximum header length check ---
                if max_header_len > 0:
                    header_content = shortened_header[1:]  # Exclude '>'
                    if len(header_content) > max_header_len:
                        log_stats["headers_exceeding_max_length"] += 1
                        log_stats["has_issues"] = True

                        # Only log details for a limited number of headers
                        if (
                            max_detailed_changes > 0
                            and len(log_stats["detailed_changes"])
                            < max_detailed_changes
                        ):
                            log_stats["detailed_changes"].append(
                                {
                                    "header": header_content,
                                    "change_type": "header_too_long",
                                    "length": len(header_content),
                                    "max_allowed": max_header_len,
                                    "details": f"Header exceeds maximum allowed length of {max_header_len} characters",
                                }
                            )

                            # Print the first few offending headers to stderr
                            if log_stats["headers_exceeding_max_length"] <= 5:
                                sys.stderr.write(
                                    f"ERROR: FASTA header exceeds maximum allowed length of {max_header_len} characters: "
                                    f"'{header_content}' (length {len(header_content)})\n"
                                )

                # Check if header was shortened
                if shortened_header != original_header:
                    log_stats["headers_shortened"] += 1
                    log_stats["has_issues"] = True
                    if (
                        max_detailed_changes > 0
                        and len(log_stats["detailed_changes"]) < max_detailed_changes
                    ):
                        log_stats["detailed_changes"].append(
                            {
                                "header": original_header[1:],
                                "change_type": "header_shortened",
                                "original": original_header[1:],
                                "modified": shortened_header[1:],
                            }
                        )

                # Check for problematic characters in the header
                has_problematic_chars = bool(re.search(r"[,;\s|:]", shortened_header))
                if has_problematic_chars:
                    log_stats["headers_with_problematic_chars"] += 1
                    log_stats["has_issues"] = True

                # Check for commas specifically
                if "," in original_header:
                    headers_with_commas += 1

                # Sanitise the header
                current_header = sanitise_header(shortened_header)

                # Log the sanitised header if it differs from the shortened header
                if current_header[1:] != shortened_header[1:]:
                    if (
                        max_detailed_changes > 0
                        and len(log_stats["detailed_changes"]) < max_detailed_changes
                    ):
                        log_stats["detailed_changes"].append(
                            {
                                "header": shortened_header[1:],
                                "change_type": "header_sanitised",
                                "original": shortened_header[1:],
                                "modified": current_header[1:],
                            }
                        )

                # Check for duplicate headers
                if current_header in headers_list:
                    if allow_duplicate_headers is False:
                        log_stats["duplicate_headers_detected"] += 1
                        log_stats["has_issues"] = True
                        if (
                            max_detailed_changes > 0
                            and len(log_stats["detailed_changes"])
                            < max_detailed_changes
                        ):
                            log_stats["detailed_changes"].append(
                                {
                                    "header": current_header[1:],
                                    "change_type": "duplicate_header",
                                    "details": f"Duplicate header found in {fasta_path}",
                                }
                            )
                        sys.stderr.write(
                            "Duplicate FASTA headers ({}) were found in the input file ({}) after truncating the headers with a delimiter\n".format(
                                current_header[1:], fasta_path
                            )
                        )
                        # Update log file before exiting if specified
                        if log_file:
                            with open(log_file, "w") as f:
                                json.dump(log_stats, f, indent=2)
                        # Use exit code 125 for validation failures - this will be used by Nextflow to avoid retries
                        sys.exit(125)
                    else:
                        log_stats["duplicate_headers_detected"] += 1
                        log_stats["has_issues"] = True
                        if (
                            max_detailed_changes > 0
                            and len(log_stats["detailed_changes"])
                            < max_detailed_changes
                        ):
                            log_stats["detailed_changes"].append(
                                {
                                    "header": current_header[1:],
                                    "change_type": "duplicate_header_allowed",
                                    "details": "Duplicate header allowed by --allow_duplicate_headers flag",
                                }
                            )

                headers_list.append(current_header)
                original_headers[current_header] = original_header
                current_sequence = []
            else:
                # Check for non-ATGC bases in the sequence
                line_upper = line.upper()
                non_atgc_count = sum(1 for c in line_upper if c not in "ATGCN\n\r")

                if non_atgc_count > 0:
                    if current_header not in log_stats.get(
                        "sequences_with_changes", {}
                    ):
                        log_stats["sequences_with_non_atgc"] += 1
                        log_stats["has_issues"] = True

                    log_stats["non_atgc_bases_replaced"] += non_atgc_count

                    # Only log the first occurrence of changes for each sequence to avoid excessive logging
                    if (
                        max_detailed_changes > 0
                        and len(log_stats["detailed_changes"]) < max_detailed_changes
                        and not any(
                            change["header"] == current_header[1:]
                            and change["change_type"] == "non_atgc_bases_replaced"
                            for change in log_stats["detailed_changes"]
                        )
                    ):
                        log_stats["detailed_changes"].append(
                            {
                                "header": current_header[1:],
                                "change_type": "non_atgc_bases_replaced",
                                "count": non_atgc_count,
                                "details": f"Non-ATGC bases replaced with N",
                            }
                        )

                # Add sanitised sequence line
                current_sequence.append(sanitise_sequence(line))

        # Process the last sequence
        process_sequence()

        # Update log with headers with commas
        if headers_with_commas > 0:
            log_stats["headers_with_commas"] = headers_with_commas
            log_stats["has_issues"] = True
            sys.stderr.write(
                "Warning: {} FASTA header(s) contained commas that were replaced with underscores\n".format(
                    headers_with_commas
                )
            )

        # Check if any headers exceeded the maximum length and exit with error if so
        if max_header_len > 0 and log_stats["headers_exceeding_max_length"] > 0:
            error_message = f"ERROR: {log_stats['headers_exceeding_max_length']} FASTA header(s) exceed maximum allowed length of {max_header_len} characters (required for FCS-adaptor)."
            sys.stderr.write(error_message + "\n")
            log_stats["error"] = error_message

            # Update log file before exiting if specified
            if log_file:
                with open(log_file, "w") as f:
                    json.dump(log_stats, f, indent=2)
            # Use exit code 125 for validation failures - this will be used by Nextflow to avoid retries
            sys.exit(125)

    # Add summary if we limited the detailed changes
    if (
        max_detailed_changes > 0
        and len(log_stats["detailed_changes"]) >= max_detailed_changes
    ):
        log_stats["detailed_changes_truncated"] = True
        log_stats["detailed_changes"].append(
            {
                "header": "SUMMARY",
                "change_type": "summary",
                "details": f"Only showing {max_detailed_changes} of the changes. Set --max_detailed_changes to a higher value to see more, or 0 to disable detailed reporting.",
            }
        )
    elif max_detailed_changes == 0:
        log_stats["detailed_changes"] = []
        log_stats["detailed_changes_truncated"] = True

    # Write log file if specified
    if log_file:
        with open(log_file, "w") as f:
            json.dump(log_stats, f, indent=2)


if __name__ == "__main__":
    args = parse_args()
    main(
        args.fasta_path,
        args.delimiter,
        args.allow_duplicate_headers,
        args.keep_n_sequences,
        args.log_file,
        args.max_detailed_changes,
        args.max_header_len,
    )
