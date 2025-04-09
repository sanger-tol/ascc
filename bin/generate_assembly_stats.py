#!/usr/bin/env python3

"""
Script to generate assembly statistics and information files for the ASCC HTML report.

This script takes an assembly FASTA file and generates:
- assembly_name.txt: Contains the name of the assembly file
- assembly_stats.txt: Contains assembly statistics (sequence count, total length, N50, GC content, etc.)
- sanitisation_report.txt: Contains information about any sanitization performed
"""

import os
import sys
import argparse
from pathlib import Path
from collections import OrderedDict
import textwrap

VERSION = "1.0.0"


def parse_args():
    parser = argparse.ArgumentParser(
        prog="ASCC Assembly Stats Generator",
        description=__doc__,
    )
    parser.add_argument("--assembly", type=str, required=True, help="Path to assembly FASTA file")
    parser.add_argument("--output_dir", type=str, default=".", help="Directory to write output files")
    parser.add_argument("--sanitisation_report", type=str, help="Path to sanitisation report file (if available)")
    parser.add_argument("-v", "--version", action="version", version=VERSION)
    return parser.parse_args()


def calculate_n50(contig_lengths):
    """Calculate N50 of a set of contig lengths."""
    if not contig_lengths:
        return 0

    sorted_lengths = sorted(contig_lengths, reverse=True)
    total_length = sum(sorted_lengths)
    target_length = total_length / 2

    current_length = 0
    for length in sorted_lengths:
        current_length += length
        if current_length >= target_length:
            return length

    return 0


def calculate_gc_content(sequence):
    """Calculate GC content of a sequence."""
    sequence = sequence.upper()
    g_count = sequence.count("G")
    c_count = sequence.count("C")
    gc_count = g_count + c_count
    total_bases = len(sequence) - sequence.count("N")

    if total_bases == 0:
        return 0

    return (gc_count / total_bases) * 100


def parse_fasta(fasta_file):
    """Parse a FASTA file and return sequence information."""
    sequences = {}
    current_seq = None
    current_header = None

    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_seq and current_header:
                    sequences[current_header] = current_seq

                current_header = line[1:].split()[0]  # Get the first word after '>'
                current_seq = ""
            else:
                current_seq += line

        # Add the last sequence
        if current_seq and current_header:
            sequences[current_header] = current_seq

    return sequences


def generate_assembly_stats(fasta_file, output_dir):
    """Generate assembly statistics and write to files."""
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Parse FASTA file
    sequences = parse_fasta(fasta_file)

    # Calculate statistics
    stats = OrderedDict()
    stats["Assembly name"] = os.path.basename(fasta_file)
    stats["Number of sequences"] = len(sequences)

    contig_lengths = [len(seq) for seq in sequences.values()]
    stats["Total length (bp)"] = sum(contig_lengths)
    stats["Longest sequence (bp)"] = max(contig_lengths) if contig_lengths else 0
    stats["Shortest sequence (bp)"] = min(contig_lengths) if contig_lengths else 0
    stats["N50 (bp)"] = calculate_n50(contig_lengths)

    # Calculate GC content
    all_sequences = "".join(sequences.values())
    stats["GC content (%)"] = f"{calculate_gc_content(all_sequences):.2f}"

    # Count Ns
    n_count = all_sequences.upper().count("N")
    stats["N count"] = n_count
    stats["N percentage (%)"] = f"{(n_count / len(all_sequences) * 100):.2f}" if all_sequences else "0.00"

    # Write assembly_name.txt
    with open(os.path.join(output_dir, "assembly_name.txt"), "w") as f:
        f.write(stats["Assembly name"])

    # Write assembly_stats.txt
    with open(os.path.join(output_dir, "assembly_stats.txt"), "w") as f:
        for key, value in stats.items():
            f.write(f"{key}: {value}\n")

    return stats


def generate_sanitisation_report(sanitisation_report_file, output_dir):
    """Generate sanitisation report or create a default one if not available."""
    os.makedirs(output_dir, exist_ok=True)

    if sanitisation_report_file and os.path.isfile(sanitisation_report_file):
        # Copy the existing sanitisation report
        with open(sanitisation_report_file, "r") as src, open(
            os.path.join(output_dir, "sanitisation_report.txt"), "w"
        ) as dst:
            dst.write(src.read())
    else:
        # Create a default sanitisation report
        with open(os.path.join(output_dir, "sanitisation_report.txt"), "w") as f:
            f.write("No sanitisation report available.\n")
            f.write("The assembly was processed by the ASCC pipeline without specific sanitisation steps.\n")


def main():
    args = parse_args()

    # Check if assembly file exists
    if not os.path.isfile(args.assembly):
        sys.stderr.write(f"Error: Assembly file '{args.assembly}' not found\n")
        sys.exit(1)

    # Generate assembly statistics
    sys.stderr.write(f"Generating assembly statistics for '{args.assembly}'...\n")
    stats = generate_assembly_stats(args.assembly, args.output_dir)

    # Generate sanitisation report
    sys.stderr.write("Generating sanitisation report...\n")
    generate_sanitisation_report(args.sanitisation_report, args.output_dir)

    sys.stderr.write(f"Assembly information files generated in '{args.output_dir}'\n")


if __name__ == "__main__":
    main()
