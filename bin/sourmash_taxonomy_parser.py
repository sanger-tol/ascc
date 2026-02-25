#!/usr/bin/env python

import sys
import argparse
import logging
import re
from collections import defaultdict
from typing import Dict
import pandas as pd
import os

# Version of sourmash_taxonomy_parser
__version__ = "1.1.0"

def setup_logger():
    """Configure and return a logger for the application."""
    logger = logging.getLogger('sourmash_taxonomy_parser')
    logger.setLevel(logging.INFO)

    # Create console handler
    ch = logging.StreamHandler(sys.stderr)
    ch.setLevel(logging.INFO)

    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)

    # Add handler to logger
    logger.addHandler(ch)
    return logger

def log_run_parameters(logger, sourmash_files, assembly_dbs, target_taxa, output_file, description):
    """Log run parameters to the logger."""
    logger.info("=" * 60)
    logger.info(f"{description} generation parameters:")
    sourmash_files_str = ', '.join(sourmash_files) if sourmash_files else 'N/A'
    logger.info(f"  Sourmash results files: {sourmash_files_str}")
    assembly_dbs_str = ', '.join(assembly_dbs) if assembly_dbs else 'N/A'
    logger.info(f"  Assembly database files: {assembly_dbs_str}")
    if target_taxa:
        target_taxa_str = ', '.join([f"{k}:{v}" for k, v in target_taxa.items()])
        logger.info(f"  Target taxa: {target_taxa_str}")
    else:
        logger.info("  Target taxa: None")
    logger.info(f"  Output file: {output_file}")
    logger.info("=" * 60)

def parse_target_taxa(target_taxa_args, logger):
    """Parse target taxa from command line arguments."""
    target_taxa = {}
    if target_taxa_args:
        for item in target_taxa_args:
            try:
                level, value = item.split(':', 1)  # Split only on first colon
                target_taxa[level.lower().strip()] = value.lower().strip()
                logger.info(f"Added target taxa filter: {level.lower().strip()} = {value.lower().strip()}")
            except ValueError:
                logger.warning(f"Invalid target taxa format '{item}'. Expected format is 'taxon:value'")
    return target_taxa

def validate_file_paths(file_paths, file_type):
    """Validate that all file paths exist."""
    for file_path in file_paths:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"{file_type} file not found: {file_path}")

def extract_contig_number(query_name):
    """
    Extract numeric contig number from query_name.

    Examples:
        'contig_1' -> 1
        'scaffold_42' -> 42
        'chr1' -> 1
        'some_name_123' -> 123
        'no_number' -> float('inf')  # sort to end if no number found

    Returns:
        int: The extracted contig number, or infinity if not found
    """
    # Look for one or more digits at the end of the string, possibly preceded by underscore or other separator
    match = re.search(r'(\d+)$', query_name)
    if match:
        return int(match.group(1))
    else:
        # No number found - sort these to the end
        return float('inf')

def parse_sourmash_results(file_paths):
    """Parse the sourmash results files and extract query-match relationships."""
    validate_file_paths(file_paths, "Sourmash results")

    query_matches = defaultdict(list)
    all_dfs = []

    logger.info(f"Processing {len(file_paths)} sourmash results file(s)")

    for file_path in file_paths:
        try:
            df = pd.read_csv(file_path)
            if df.empty:
                logger.warning(f"Sourmash results file {file_path} is empty")
                continue
            required_columns = ['query_name', 'match_name', 'containment', 'jaccard']
            missing_columns = [col for col in required_columns if col not in df.columns]
            if missing_columns:
                raise ValueError(f"Missing required columns in {file_path}: {missing_columns}")
            logger.info(f"Loaded {len(df)} sourmash results from {file_path}")
            all_dfs.append(df)
        except pd.errors.EmptyDataError:
            logger.warning(f"Sourmash results file {file_path} is empty or has no columns")
            continue
        except Exception as e:
            logger.error(f"Error loading sourmash results file {file_path}: {e}")
            raise

    if not all_dfs:
        raise ValueError("No valid sourmash results files provided")

    # Concatenate all DataFrames
    combined_df = pd.concat(all_dfs, ignore_index=True)
    logger.info(f"Combined sourmash results: {len(combined_df)} total rows")

    for _, row in combined_df.iterrows():
        query_name = row['query_name']
        match_name = row['match_name']

        # Skip rows where containment or jaccard are not numeric (duplicate headers)
        try:
            containment = float(row['containment'])
            jaccard = float(row['jaccard'])
            # Ensure intersect_hashes is numeric (float or int)
            intersect_hashes = float(row.get('intersect_hashes', 0)) if pd.notna(row.get('intersect_hashes', 0)) else 0.0
        except (ValueError, TypeError):
            # Skip this row - it's likely a duplicate header
            logger.warning(f"Skipping row with non-numeric values: query_name={query_name}, match_name={match_name}")
            continue

        if " " in match_name:
            # Handle cases where match_name contains spaces
            match_name = match_name.strip().split(" ")[0]

        query_matches[query_name].append({
            'match_name': match_name,
            'containment': containment,
            'jaccard': jaccard,
            'intersect_hashes': intersect_hashes
        })

    return query_matches

def parse_assembly_database(file_paths):
    """Parse the assembly database files and return a DataFrame with taxonomic information."""
    validate_file_paths(file_paths, "Assembly database")

    all_dfs = []

    logger.info(f"Processing {len(file_paths)} assembly database file(s)")

    for file_path in file_paths:
        try:
            assembly_df = pd.read_csv(file_path)
            if assembly_df.empty:
                logger.warning(f"Assembly database file {file_path} is empty")
                continue
            logger.info(f"Loaded {len(assembly_df)} assembly records from {file_path}")
            all_dfs.append(assembly_df)
        except pd.errors.EmptyDataError:
            logger.warning(f"Assembly database file {file_path} is empty or has no columns")
            continue
        except Exception as e:
            logger.error(f"Error loading assembly database file {file_path}: {e}")
            raise

    if not all_dfs:
        raise ValueError("No valid assembly database files provided")

    # Concatenate all DataFrames
    combined_df = pd.concat(all_dfs, ignore_index=True)
    logger.info(f"Combined assembly database: {len(combined_df)} total rows")

    # Check for duplicate indices if setting assembly_accession
    if 'assembly_accession' in combined_df.columns:
        if combined_df['assembly_accession'].duplicated().any():
            logger.warning("Duplicate assembly_accession values found after merging. Keeping first occurrence.")
            combined_df = combined_df.drop_duplicates(subset='assembly_accession')
        combined_df.set_index('assembly_accession', inplace=True)
    else:
        # Use first column as index if standard names not found
        first_col = combined_df.columns[0]
        if combined_df[first_col].duplicated().any():
            logger.warning(f"Duplicate values in first column '{first_col}' after merging. Keeping first occurrence.")
            combined_df = combined_df.drop_duplicates(subset=first_col)
        combined_df.set_index(first_col, inplace=True)

    # Filter out assemblies without valid taxid
    if 'taxid' in combined_df.columns:
        initial_count = len(combined_df)
        combined_df = combined_df[combined_df['taxid'].notna() & (combined_df['taxid'].astype(str).str.strip() != 'nan')]
        filtered_count = len(combined_df)
        logger.info(f"Filtered out {initial_count - filtered_count} assemblies without valid taxid")

    return combined_df

def generate_summary(query_matches):
    """Generate summary of taxonomic matches for each query."""
    summary = []

    for query_name, matches in query_matches.items():
        # Sort matches by intersect_hashes in descending order (highest first)
        sorted_matches = sorted(matches, key=lambda x: x['intersect_hashes'], reverse=True)

        # Generate rows with rank
        for rank, match in enumerate(sorted_matches, 1):
            summary.append((
                query_name,
                match['match_name'],
                rank,
                match['containment'],
                match['jaccard'],
                match['intersect_hashes']
            ))

    return summary

def get_target_genomes(assembly_df, target_taxa):
    """Pre-compute set of all genomes that match target taxa."""
    if assembly_df is None or not target_taxa:
        return set()

    target_genomes = set()

    for level, target_value in target_taxa.items():
        if level in assembly_df.columns:
            # Vectorized search for this taxonomic level
            mask = (
                assembly_df[level].notna() &
                (assembly_df[level].astype(str).str.lower().str.strip() == target_value) &
                (assembly_df[level].astype(str).str.strip() != 'nan')
            )
            level_matches = set(assembly_df[mask].index)
            target_genomes.update(level_matches)

            logger.info(f"Found {len(level_matches)} genomes matching {level}={target_value}")

    logger.info(f"Total target genomes: {len(target_genomes)}")
    return target_genomes

def write_summary_output(summary, output_file, assembly_df=None, target_taxa=None, sourmash_files=None, assembly_dbs=None):
    """Write summary output to a file."""
    # Pre-compute target genomes set
    target_genomes = get_target_genomes(assembly_df, target_taxa)

    # Pre-compute taxid dictionary for faster lookup
    taxid_dict = {}
    if assembly_df is not None and 'taxid' in assembly_df.columns:
        taxid_dict = assembly_df['taxid'].dropna().astype(str).str.strip().to_dict()

    # Sort summary with priority:
    # 1. By contig number in query_name (ascending: contig_1, contig_2, ...)
    # 2. By intersect_hashes (descending: more hashes first)
    sorted_summary = sorted(
        summary,
        key=lambda x: (extract_contig_number(x[0]), -float(x[5]))
    )

    # Log run parameters
    log_run_parameters(logger, sourmash_files, assembly_dbs, target_taxa, output_file, "Summary file")

    with open(output_file, 'w') as f:
        # Write header only (no comments)
        f.write("header,assembly_accession,taxa,top_n,containment,jaccard,intersect_hashes,is_target\n")

        for row in sorted_summary:
            query_name, match_name, rank, containment, jaccard, intersect_hashes = row

            # Fast lookup: is this genome in target set?
            is_target = match_name in target_genomes if target_genomes else None

            # Fast lookup: get taxid from pre-computed dictionary
            taxa = taxid_dict.get(match_name, "nan")
            if taxa == "nan":
                logger.warning(f"Assembly {match_name} not found in assembly_taxa database")

            # Format the values
            is_target_str = str(is_target) if is_target is not None else "None"

            # Write the output
            f.write(f"{query_name},{match_name},{taxa},{rank},{containment:.6f},{jaccard:.6f},{intersect_hashes},{is_target_str}\n")

def write_non_target_output(summary_file, output_file, assembly_df, sourmash_files=None, assembly_dbs=None, target_taxa=None):
    """Write non-target queries with lineage information of top1 match."""
    if not os.path.exists(summary_file):
        raise FileNotFoundError(f"Summary file not found: {summary_file}")

    # Read the summary file to identify non-target queries (no comment lines now)
    df = pd.read_csv(summary_file)

    # Group by header and check if any row has is_target=True
    query_groups = df.groupby('header')
    non_target_queries = []

    for query_name, group in query_groups:
        # Check if this query has any target taxa matches
        has_target = group['is_target'].any()

        if not has_target:
            # Get the top1 match (top_n == 1)
            top1_rows = group[group['top_n'] == 1]
            if not top1_rows.empty:
                top1_row = top1_rows.iloc[0]
                taxa = top1_row['taxa']
                # Skip if taxa is None/nan (unknown taxonomy)
                if pd.isna(taxa) or str(taxa).strip() in ['None', 'nan']:
                    continue
                non_target_queries.append(top1_row)

    # Log run parameters
    log_run_parameters(logger, sourmash_files, assembly_dbs, target_taxa, output_file, "Non-target file")
    logger.info(f"  Non-target queries found: {len(non_target_queries)}")

    # Write the non-target queries file
    with open(output_file, 'w') as f:
        # Write header only (no comments)
        f.write("header,assembly_accession,taxa,containment,jaccard,intersect_hashes,species,genus,family,order,class,phylum,kingdom\n")

        for _, row in enumerate(non_target_queries):
            query_name = row['header']
            match_name = row['assembly_accession']
            taxa = row['taxa']
            containment = row['containment']
            jaccard = row['jaccard']
            intersect_hashes = row['intersect_hashes']

            # Get lineage information using vectorized operations
            lineage_info = ["None"] * 7  # species, genus, family, order, class, phylum, kingdom
            try:
                if match_name in assembly_df.index:
                    lineage_row = assembly_df.loc[match_name]
                    lineage_columns = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']
                    for i, col in enumerate(lineage_columns):
                        if col in lineage_row.index and pd.notna(lineage_row[col]) and str(lineage_row[col]).strip() != 'nan':
                            lineage_info[i] = str(lineage_row[col]).strip()
            except Exception as e:
                logger.debug(f"Error getting lineage for {match_name}: {e}")

            # Write the output
            lineage_str = ",".join(lineage_info)
            f.write(f"{query_name},{match_name},{taxa},{containment:.6f},{jaccard:.6f},{intersect_hashes},{lineage_str}\n")

def main(sourmash_files=None, assembly_dbs=None, target_taxa=None, outdir=None):
    """Main function to run the sourmash parser."""

    # Initialize variables
    assembly_df = None
    query_matches = defaultdict(list)

    # Process assembly database if provided
    if assembly_dbs:
        logger.info(f"Processing assembly database files: {assembly_dbs}")
        assembly_df = parse_assembly_database(assembly_dbs)

    logger.info(f"Processing sourmash results files: {sourmash_files}")
    query_matches = parse_sourmash_results(sourmash_files)

    # Generate summary
    summary = generate_summary(query_matches)

    # Create output directory if it doesn't exist
    os.makedirs(outdir, exist_ok=True)

    # Generate output filenames
    # Use 'sourmash_results' as basename if multiple files, else use the single file's basename
    sourmash_basename = "sourmash_results" if len(sourmash_files) > 1 else os.path.splitext(os.path.basename(sourmash_files[0]))[0]

    # Create target taxa suffix for filename
    taxa_suffix = ""
    if target_taxa:
        taxa_parts = []
        for level, value in target_taxa.items():
            taxa_parts.append(f"{level}_{value}")
        taxa_suffix = "_" + "_".join(taxa_parts)

    summary_filename = f"{sourmash_basename}{taxa_suffix}.summary.csv"
    summary_file = os.path.join(outdir, summary_filename)

    logger.info(f"Writing summary to: {summary_file}")

    # Write summary output
    write_summary_output(summary, summary_file, assembly_df, target_taxa, sourmash_files, assembly_dbs)

    # Write non-target queries if assembly database and target taxa are available
    if assembly_df is not None and target_taxa:
        non_target_filename = f"{sourmash_basename}{taxa_suffix}.non_target.csv"
        non_target_file = os.path.join(outdir, non_target_filename)
        logger.info(f"Writing non-target queries to: {non_target_file}")
        write_non_target_output(summary_file, non_target_file, assembly_df, sourmash_files, assembly_dbs, target_taxa)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse sourmash results and extract query-match relationships with similarity scores.")
    parser.add_argument('-s', '--sourmash_results', nargs='+', help='Path(s) to the sourmash results file(s) to parse', required=False)
    parser.add_argument('-a', '--assembly_db', nargs='+', help='Path(s) to assembly database file(s) with taxonomic information', required=False)
    parser.add_argument('--target_taxa', nargs='+', help='Target taxa in format taxon:value (e.g., order:coleoptera family:Carabidae)')
    parser.add_argument('--log', help='Path to log file (if not provided, logs to stderr)', default=None)
    parser.add_argument('-o', '--outdir', help='Output directory for results (default: current directory)', default=None)
    parser.add_argument('--version', action='version', version=f'sourmash_taxonomy_parser {__version__}')
    args = parser.parse_args()

    # If only version was requested, exit
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    # Validate required arguments (unless --version was used)
    if not args.sourmash_results or not args.assembly_db or not args.outdir:
        parser.error("the following arguments are required: -s/--sourmash_results, -a/--assembly_db, -o/--outdir")

    # Initialize logger
    logger = setup_logger()

    # Configure file logging if requested
    if args.log:
        file_handler = logging.FileHandler(args.log)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)

    # Parse target taxa
    target_taxa = parse_target_taxa(args.target_taxa, logger)

    main(sourmash_files=args.sourmash_results, assembly_dbs=args.assembly_db, target_taxa=target_taxa, outdir=args.outdir)
