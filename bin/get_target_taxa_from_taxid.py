#!/usr/bin/env python3

"""
Script to extract target taxa from a taxid using NCBI rankedlineage.dmp file.
Returns target_taxa in format 'level:taxa_name' (e.g., 'order:Artiodactyla').

By default, writes the full fallback chain starting from the requested taxonomy
level up to kingdom, one entry per line.  The downstream parser
(sourmash_taxonomy_parser.py) will then pick the first level that has enough
genomes in the assembly database (see --min_target_genomes).

Example output for --level order when the full lineage is available:
    order:Artiodactyla
    class:Mammalia
    phylum:Chordata
    kingdom:Metazoa

Usage:
    get_target_taxa_from_taxid.py --taxid <taxid> --rankedlineage <path> --level <level> --output <output_file>

Author: Danil Zilov (@zilov)
Version: 1.1.0
"""

import argparse
import sys
import os


VERSION = "1.1.0"

# Mapping of taxonomy levels to their column indices in rankedlineage.dmp
# Format: taxid | name | species | genus | family | order | class | phylum | kingdom | domain
TAXONOMY_LEVELS = {
    'species': 2,
    'genus': 3,
    'family': 4,
    'order': 5,
    'class': 6,
    'phylum': 7,
    'kingdom': 8,
    'domain': 9
}

# Ordered fallback chain: from most specific to least specific.
# Used to build the fallback list written to the output file.
FALLBACK_ORDER = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']


def parse_rankedlineage(rankedlineage_path, target_taxid):
    """
    Parse the rankedlineage.dmp file and find the lineage for the given taxid.

    Args:
        rankedlineage_path: Path to the rankedlineage.dmp file
        target_taxid: Target taxid to search for

    Returns:
        Dictionary with taxonomy information or None if not found
    """
    if not os.path.isfile(rankedlineage_path):
        sys.stderr.write(
            f"ERROR: The NCBI rankedlineage.dmp file was not found at: {rankedlineage_path}\n"
        )
        sys.exit(1)

    with open(rankedlineage_path, 'r') as f:
        for line in f:
            split_line = line.strip().split('|')
            split_line = [n.strip() for n in split_line]

            # Validate line format
            if len(split_line) < 10:
                continue

            taxid = split_line[0]
            if taxid == target_taxid:
                return {
                    'taxid': split_line[0],
                    'name': split_line[1],
                    'species': split_line[2],
                    'genus': split_line[3],
                    'family': split_line[4],
                    'order': split_line[5],
                    'class': split_line[6],
                    'phylum': split_line[7],
                    'kingdom': split_line[8],
                    'domain': split_line[9]
                }

    return None


def get_target_taxa(lineage_dict, taxonomy_level):
    """
    Build a fallback chain of target taxa starting from taxonomy_level up to kingdom.

    Args:
        lineage_dict: Dictionary with taxonomy information
        taxonomy_level: Taxonomy level to start from (e.g., 'order', 'class')

    Returns:
        List of strings in format 'level:taxa_name', from most specific (taxonomy_level)
        to least specific (kingdom).  Empty strings / missing levels are skipped.
        Returns an empty list if the requested level itself is empty.
    """
    if taxonomy_level not in TAXONOMY_LEVELS:
        sys.stderr.write(
            f"ERROR: Invalid taxonomy level '{taxonomy_level}'. "
            f"Valid levels are: {', '.join(TAXONOMY_LEVELS.keys())}\n"
        )
        sys.exit(1)

    # Start from the requested level; include all levels above it up to kingdom
    try:
        start_idx = FALLBACK_ORDER.index(taxonomy_level)
    except ValueError:
        # taxonomy_level not in FALLBACK_ORDER (e.g. 'domain') — just use that single level
        start_idx = len(FALLBACK_ORDER)

    chain = []
    # Always try the requested level first; if it's empty, return empty list (caller handles)
    primary_value = lineage_dict.get(taxonomy_level, '').strip()
    if not primary_value:
        return []

    for level in FALLBACK_ORDER[start_idx:]:
        value = lineage_dict.get(level, '').strip()
        if value:
            chain.append(f"{level}:{value}")

    return chain


def parse_args():
    """
    Parse command line arguments.

    Returns:
        Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Extract target taxa from taxid using NCBI rankedlineage.dmp",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '--taxid',
        required=True,
        help='Target taxid to search for'
    )

    parser.add_argument(
        '--rankedlineage',
        required=True,
        help='Path to NCBI rankedlineage.dmp file'
    )

    parser.add_argument(
        '--level',
        required=True,
        choices=list(TAXONOMY_LEVELS.keys()),
        help='Taxonomy level to extract (e.g., order, class, family)'
    )

    parser.add_argument(
        '--output',
        required=True,
        help='Output file path for target_taxa'
    )

    parser.add_argument(
        '--version',
        action='version',
        version=f'%(prog)s {VERSION}'
    )

    return parser.parse_args()


def main():
    args = parse_args()

    # Parse rankedlineage file
    lineage_dict = parse_rankedlineage(args.rankedlineage, args.taxid)

    if lineage_dict is None:
        sys.stderr.write(
            f"WARNING: Taxid '{args.taxid}' not found in rankedlineage file.\n"
            "Sourmash will be skipped for this sample.\n"
        )
        # Write empty file to indicate failure but allow pipeline to continue
        with open(args.output, 'w') as f:
            pass  # Create empty file
        sys.exit(0)

    # Build fallback chain from requested level upwards
    target_taxa_chain = get_target_taxa(lineage_dict, args.level)

    if not target_taxa_chain:
        sys.stderr.write(
            f"WARNING: Taxonomy level '{args.level}' is empty for taxid '{args.taxid}'.\n"
            f"Available lineage: {lineage_dict}\n"
            f"Sourmash will be skipped for this sample.\n"
        )
        # Write empty file to indicate failure but allow pipeline to continue
        with open(args.output, 'w') as f:
            pass  # Create empty file
        sys.exit(0)

    # Write fallback chain — comma-separated on a single line, most specific first
    with open(args.output, 'w') as f:
        f.write(','.join(target_taxa_chain) + '\n')

    sys.stdout.write(
        f"Successfully extracted target taxa chain ({len(target_taxa_chain)} levels): "
        f"{', '.join(target_taxa_chain)}\n"
    )


if __name__ == '__main__':
    main()
