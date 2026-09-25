#!/usr/bin/env python3

"""
Script to extract target taxa from a taxid using NCBI rankedlineage.dmp file.
Returns target_taxa in format 'level:taxa_name' (e.g., 'order:Artiodactyla').
With --fallback_levels, the higher levels follow on the same line, space-separated, in order
(e.g. 'order:Andreaeales class:Andreaeopsida phylum:Streptophyta kingdom:Viridiplantae'); empty
levels are skipped. sourmash_taxonomy_parser.py --target_fallback uses the first one present in its DB.

Usage:
    get_target_taxa_from_taxid.py --taxid <taxid> --rankedlineage <path> --level <level> --output <output_file>
                                  [--fallback_levels class,phylum,kingdom]

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
    Extract target taxa from lineage dictionary at specified taxonomy level.

    Args:
        lineage_dict: Dictionary with taxonomy information
        taxonomy_level: Taxonomy level to extract (e.g., 'order', 'class')

    Returns:
        String in format 'level:taxa_name' or None if level is empty
    """
    if taxonomy_level not in TAXONOMY_LEVELS:
        sys.stderr.write(
            f"ERROR: Invalid taxonomy level '{taxonomy_level}'. "
            f"Valid levels are: {', '.join(TAXONOMY_LEVELS.keys())}\n"
        )
        sys.exit(1)

    taxa_name = lineage_dict.get(taxonomy_level, '')

    # Check if the taxa name is empty or just whitespace
    if not taxa_name or taxa_name.strip() == '':
        return None

    return f"{taxonomy_level}:{taxa_name}"


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
        '--fallback_levels',
        default='',
        help='Comma-separated higher levels to append, in order (e.g. class,phylum,kingdom)'
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

    # Extract target taxa at specified level
    target_taxa = get_target_taxa(lineage_dict, args.level)
    # only levels above the requested one (e.g. --level phylum ignores class)
    fallback_levels = [lv.strip() for lv in args.fallback_levels.split(',')
                       if lv.strip() in TAXONOMY_LEVELS and TAXONOMY_LEVELS[lv.strip()] > TAXONOMY_LEVELS[args.level]]
    if target_taxa is None and fallback_levels:
        # e.g. no order assigned in NCBI: start from the first non-empty fallback level
        for level in fallback_levels:
            target_taxa = get_target_taxa(lineage_dict, level)
            if target_taxa:
                sys.stderr.write(f"WARNING: level '{args.level}' is empty for taxid '{args.taxid}'; starting at {target_taxa}\n")
                break

    if target_taxa is None:
        sys.stderr.write(
            f"WARNING: Taxonomy level '{args.level}' is empty for taxid '{args.taxid}'.\n"
            f"Available lineage: {lineage_dict}\n"
            f"Sourmash will be skipped for this sample.\n"
        )
        # Write empty file to indicate failure but allow pipeline to continue
        with open(args.output, 'w') as f:
            pass  # Create empty file
        sys.exit(0)

    # Append higher levels for the parser's fallback (skipping empty or repeated levels)
    taxa_list = [target_taxa]
    for level in fallback_levels:
        extra = get_target_taxa(lineage_dict, level)
        if extra and extra not in taxa_list:
            taxa_list.append(extra)
    target_taxa = " ".join(taxa_list)

    # Write successful result
    with open(args.output, 'w') as f:
        f.write(f"{target_taxa}\n")

    sys.stdout.write(f"Successfully extracted target taxa: {target_taxa}\n")


if __name__ == '__main__':
    main()
