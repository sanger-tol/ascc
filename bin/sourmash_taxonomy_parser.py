#!/usr/bin/env python

import sys
import argparse
import logging
import re
from collections import defaultdict
from typing import Dict
import polars as pl
import os

# Version of sourmash_taxonomy_parser
__version__ = "1.3.0"

# Module-level logger — always available when the module is imported.
# Handlers are added by setup_logger() at runtime (CLI) or by the caller (tests).
logger = logging.getLogger(__name__)


def setup_logger(log_file=None):
    """Add stderr (and optionally file) handlers to the module logger."""
    logger.setLevel(logging.INFO)

    # Avoid adding duplicate handlers if called more than once
    if not logger.handlers:
        ch = logging.StreamHandler(sys.stderr)
        ch.setLevel(logging.INFO)
        ch.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
        logger.addHandler(ch)

    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
        logger.addHandler(fh)

def log_run_parameters(sourmash_files, assembly_dbs, target_taxa, output_file, description):
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

def parse_target_taxa(target_taxa_args):
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

def parse_sourmash_results(file_paths, exclude_accessions=None):
    """Parse the sourmash results files and extract query-match relationships."""
    validate_file_paths(file_paths, "Sourmash results")

    query_matches = defaultdict(list)
    all_dfs = []

    logger.info(f"Processing {len(file_paths)} sourmash results file(s)")

    # Prepare exclude set for fast lookup
    exclude_set = set()
    if exclude_accessions:
        for acc_arg in exclude_accessions:
            # Split by comma in case comma-separated values are passed
            for acc in acc_arg.split(','):
                exclude_set.add(acc.strip())
        logger.info(f"Excluding {len(exclude_set)} accessions: {', '.join(sorted(exclude_set))}")

    for file_path in file_paths:
        try:
            df = pl.read_csv(file_path, infer_schema_length=0)  # read all as Utf8 to handle duplicate headers
            if df.is_empty():
                logger.warning(f"Sourmash results file {file_path} is empty")
                continue
            required_columns = ['query_name', 'match_name', 'containment', 'jaccard']
            missing_columns = [col for col in required_columns if col not in df.columns]
            if missing_columns:
                raise ValueError(f"Missing required columns in {file_path}: {missing_columns}")
            logger.info(f"Loaded {len(df)} sourmash results from {file_path}")
            all_dfs.append(df)
        except pl.exceptions.NoDataError:
            logger.warning(f"Sourmash results file {file_path} is empty or has no columns")
            continue
        except Exception as e:
            logger.error(f"Error loading sourmash results file {file_path}: {e}")
            raise

    if not all_dfs:
        raise ValueError("No valid sourmash results files provided")

    # Concatenate all DataFrames
    combined_df = pl.concat(all_dfs, how="diagonal_relaxed")
    logger.info(f"Combined sourmash results: {len(combined_df)} total rows")

    # Drop rows where containment/jaccard are non-numeric (duplicate header rows)
    combined_df = combined_df.filter(
        pl.col("containment").str.contains(r"^[0-9]") &
        pl.col("jaccard").str.contains(r"^[0-9]")
    )

    # Cast numeric columns
    combined_df = combined_df.with_columns([
        pl.col("containment").cast(pl.Float64),
        pl.col("jaccard").cast(pl.Float64),
        pl.when(pl.col("intersect_hashes").is_not_null())
          .then(pl.col("intersect_hashes").cast(pl.Float64))
          .otherwise(pl.lit(0.0))
          .alias("intersect_hashes")
        if "intersect_hashes" in combined_df.columns
        else pl.lit(0.0).alias("intersect_hashes")
    ])

    # Strip accession from match_name (take first token before space)
    combined_df = combined_df.with_columns(
        pl.col("match_name").str.split(" ").list.first().alias("match_name")
    )

    # Exclude accessions
    filtered_count = 0
    if exclude_set:
        before = len(combined_df)
        combined_df = combined_df.filter(~pl.col("match_name").is_in(list(exclude_set)))
        filtered_count = before - len(combined_df)

    # Build query_matches dict from polars rows (avoid iterrows overhead)
    for row in combined_df.iter_rows(named=True):
        query_matches[row["query_name"]].append({
            "match_name": row["match_name"],
            "containment": row["containment"],
            "jaccard": row["jaccard"],
            "intersect_hashes": row["intersect_hashes"],
        })

    # Log filtering statistics
    total_matches = sum(len(matches) for matches in query_matches.values())
    logger.info(f"Total matches after filtering: {total_matches}")
    if exclude_set:
        logger.info(f"Filtered out {filtered_count} matches due to excluded accessions")

    return query_matches

def parse_assembly_database(file_paths):
    """Parse the assembly database files and return a polars DataFrame with taxonomic information.

    Returns a polars DataFrame. The accession column is kept as a regular column named
    'assembly_accession' (or first column) — lookups use polars filter/join instead of
    pandas index.
    """
    validate_file_paths(file_paths, "Assembly database")

    all_dfs = []

    logger.info(f"Processing {len(file_paths)} assembly database file(s)")

    for file_path in file_paths:
        try:
            assembly_df = pl.read_csv(file_path, infer_schema_length=1000)
            if assembly_df.is_empty():
                logger.warning(f"Assembly database file {file_path} is empty")
                continue
            logger.info(f"Loaded {len(assembly_df)} assembly records from {file_path}")
            all_dfs.append(assembly_df)
        except pl.exceptions.NoDataError:
            logger.warning(f"Assembly database file {file_path} is empty or has no columns")
            continue
        except Exception as e:
            logger.error(f"Error loading assembly database file {file_path}: {e}")
            raise

    if not all_dfs:
        raise ValueError("No valid assembly database files provided")

    # Concatenate all DataFrames
    combined_df = pl.concat(all_dfs, how="diagonal_relaxed")
    logger.info(f"Combined assembly database: {len(combined_df)} total rows")

    # Validate required columns
    required_cols = ['assembly_accession', 'taxid']
    missing_required = [c for c in required_cols if c not in combined_df.columns]
    if missing_required:
        raise ValueError(
            f"Assembly database is missing required column(s): {missing_required}. "
            f"Found columns: {combined_df.columns}"
        )

    # Deduplicate on accession column
    before = len(combined_df)
    combined_df = combined_df.unique(subset=['assembly_accession'], keep="first")
    if len(combined_df) < before:
        logger.warning("Duplicate 'assembly_accession' values found after merging. Keeping first occurrence.")

    # Filter out assemblies without valid taxid
    initial_count = len(combined_df)
    combined_df = combined_df.filter(
        pl.col('taxid').is_not_null() &
        (pl.col('taxid').cast(pl.Utf8).str.strip_chars() != 'nan')
    )
    logger.info(f"Filtered out {initial_count - len(combined_df)} assemblies without valid taxid")

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
        if level not in assembly_df.columns:
            logger.warning(
                f"Target taxa level '{level}' not found in assembly database columns: "
                f"{assembly_df.columns} — this filter will match nothing"
            )
            continue
        # Vectorized search for this taxonomic level (polars)
        level_matches = (
            assembly_df
            .filter(
                pl.col(level).is_not_null() &
                (pl.col(level).cast(pl.Utf8).str.to_lowercase().str.strip_chars() == target_value) &
                (pl.col(level).cast(pl.Utf8).str.strip_chars() != 'nan')
            )
            .get_column('assembly_accession')
            .to_list()
        )
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
        taxid_dict = {
            row[0]: row[1]
            for row in assembly_df
            .filter(pl.col('taxid').is_not_null())
            .select(['assembly_accession', pl.col('taxid').cast(pl.Utf8).str.strip_chars()])
            .iter_rows()
        }

    # Sort summary with priority:
    # 1. By contig number in query_name (ascending: contig_1, contig_2, ...)
    # 2. By intersect_hashes (descending: more hashes first)
    sorted_summary = sorted(
        summary,
        key=lambda x: (extract_contig_number(x[0]), -float(x[5]))
    )

    # Log run parameters
    log_run_parameters(sourmash_files, assembly_dbs, target_taxa, output_file, "Summary file")

    with open(output_file, 'w') as f:
        # Write header only (no comments)
        f.write("header,assembly_accession,taxa,top_n,containment,jaccard,intersect_hashes,is_target\n")

        for row in sorted_summary:
            query_name, match_name, rank, containment, jaccard, intersect_hashes = row

            # Fast lookup: get taxid from pre-computed dictionary
            taxa = taxid_dict.get(match_name, "nan")
            unknown_taxid = (taxa == "nan")
            if unknown_taxid:
                logger.warning(
                    f"Assembly {match_name} not found in assembly_taxa database — "
                    f"treating as target (safe side) to avoid accidental removal"
                )

            # Fast lookup: is this genome in target set?
            # Unknown taxonomy is treated as target (conservative: never remove unknowns automatically)
            if unknown_taxid:
                is_target = True
            elif target_genomes:
                is_target = match_name in target_genomes
            else:
                is_target = None

            # Format the values
            is_target_str = str(is_target) if is_target is not None else "None"

            # Write the output
            f.write(f"{query_name},{match_name},{taxa},{rank},{containment:.6f},{jaccard:.6f},{intersect_hashes},{is_target_str}\n")

def write_non_target_output(summary_file, output_file, assembly_df, sourmash_files=None, assembly_dbs=None, target_taxa=None):
    """Write non-target queries with lineage information of top1 match."""
    if not os.path.exists(summary_file):
        raise FileNotFoundError(f"Summary file not found: {summary_file}")

    # Read the summary file (polars)
    df = pl.read_csv(summary_file, infer_schema_length=1000)

    # Cast is_target to boolean-compatible: 'True' -> True
    df = df.with_columns(
        pl.col('is_target').cast(pl.Utf8).str.to_lowercase().str.strip_chars().alias('is_target_str')
    )

    # For each query (header), find those where NO row has is_target=True
    # and grab the top_n==1 row
    has_target_df = (
        df
        .filter(pl.col('is_target_str') == 'true')
        .select('header')
        .unique()
    )
    has_target_set = set(has_target_df.get_column('header').to_list())

    non_target_df = (
        df
        .filter(
            ~pl.col('header').is_in(list(has_target_set)) &
            (pl.col('top_n') == 1)
        )
    )

    non_target_queries = non_target_df.to_dicts()

    # Log run parameters
    log_run_parameters(sourmash_files, assembly_dbs, target_taxa, output_file, "Non-target file")
    logger.info(f"  Non-target queries found: {len(non_target_queries)}")

    # Build lineage lookup dict: accession -> {col: value}
    lineage_columns = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']
    avail_lineage_cols = [c for c in lineage_columns if c in assembly_df.columns]
    lineage_lookup: dict[str, dict] = {}
    for row in assembly_df.select(['assembly_accession'] + avail_lineage_cols).iter_rows(named=True):
        lineage_lookup[row['assembly_accession']] = row

    # Write the non-target queries file
    with open(output_file, 'w') as f:
        # Write header only (no comments)
        f.write("header,assembly_accession,taxa,containment,jaccard,intersect_hashes,species,genus,family,order,class,phylum,kingdom\n")

        for row in non_target_queries:
            query_name = row['header']
            match_name = row['assembly_accession']
            taxa = row['taxa']
            containment = float(row['containment'])
            jaccard = float(row['jaccard'])
            intersect_hashes = row['intersect_hashes']

            lineage_info = ["None"] * 7
            try:
                if match_name in lineage_lookup:
                    lrow = lineage_lookup[match_name]
                    for i, col in enumerate(lineage_columns):
                        val = lrow.get(col)
                        if val is not None and str(val).strip() not in ('nan', 'None', ''):
                            lineage_info[i] = str(val).strip()
            except Exception as e:
                logger.debug(f"Error getting lineage for {match_name}: {e}")

            lineage_str = ",".join(lineage_info)
            f.write(f"{query_name},{match_name},{taxa},{containment:.6f},{jaccard:.6f},{intersect_hashes},{lineage_str}\n")

def main(sourmash_files=None, assembly_dbs=None, target_taxa=None, outdir=None, exclude_accessions=None):
    """Main function to run the sourmash parser."""

    # Initialize variables
    assembly_df = None
    query_matches = defaultdict(list)

    # Process assembly database if provided
    if assembly_dbs:
        logger.info(f"Processing assembly database files: {assembly_dbs}")
        assembly_df = parse_assembly_database(assembly_dbs)

    logger.info(f"Processing sourmash results files: {sourmash_files}")
    query_matches = parse_sourmash_results(sourmash_files, exclude_accessions)

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
    parser.add_argument('--accessions_to_exclude', nargs='+', help='List of assembly accessions to exclude from results (comma-separated or multiple values)', default=None)
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

    # Initialise logger (add handlers; file handler only if --log provided)
    setup_logger(log_file=args.log)

    # Parse target taxa
    target_taxa = parse_target_taxa(args.target_taxa)

    main(sourmash_files=args.sourmash_results, assembly_dbs=args.assembly_db, target_taxa=target_taxa, outdir=args.outdir, exclude_accessions=args.accessions_to_exclude)
