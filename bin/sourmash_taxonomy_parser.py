#!/usr/bin/env python
"""
Parse sourmash (branchwater) multisearch results into per-hit and per-contig tables.

Outputs (in --outdir, <base> = input basename or 'sourmash_results', <taxa> = e.g. '_order_coleoptera'):
    <base><taxa>.summary.csv         one row per (contig, DB hit), ranked per contig by containment (top_n)
    <base><taxa>.classification.csv  one row per contig (all contigs with --contig_sizes), in the layout of
                                     ascc's parsed FCS-GX table: scaff, sourmash_top_tax_name, sourmash_top_taxid,
                                     sourmash_div, sourmash_action (EXCLUDE / REVIEW / NA), lineage, evidence
    <base><taxa>.non_target.csv      contigs with action EXCLUDE and their best non-target hit (backward compat)

Decision rule (per contig), see --help for the parameters:
    * hits from excluded accessions are dropped first (--accessions_to_exclude,
      --exclude_accessions_file, --exclude_tolid, --exclude_species);
    * hits are ranked by containment (fraction of the contig's hashes shared; scale-invariant, so DBs
      with different scaled are comparable), and hash counts are normalised to scaled=200 equivalents
      (norm_hashes = intersect_hashes * scaled / 200);
    * a non-target hit is *evidence* only if norm_hashes >= --min_hashes and
      containment >= --min_containment; target hits always count (conservative);
    * EXCLUDE: best evidence non-target containment > --margin * best target containment
               (or there is no target hit at all);
    * REVIEW : best evidence non-target hit beats the best target hit, but not by the margin;
    * NA     : target wins, only weak non-target evidence (like FCS-GX 'low-coverage'), no hits.
The Aug-2025 behaviour (any target hit keeps the contig, any single non-target hash removes it) is
--min_hashes 0 --min_containment 0 --margin inf.
"""

import sys
import argparse
import logging
import math
import re
import os
import polars as pl

# Version of sourmash_taxonomy_parser
__version__ = "2.2.0"

# Module-level logger — always available when the module is imported.
# Handlers are added by setup_logger() at runtime (CLI) or by the caller (tests).
logger = logging.getLogger(__name__)

REFERENCE_SCALED = 200
LINEAGE_COLUMNS = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']
# A Tree of Life ID: lowercase prefix, Genus (3+), species (3+), specimen number, e.g. icCatNigr1, cbHetHete15
TOLID_RE = re.compile(r'^[a-z]{1,3}[A-Z][a-z]{2,}[A-Z][a-z]{2,}\d+$')

# Defaults chosen on the 18-sample F1 harness + icCatNigr1 (TODO 04 results):
# 1 shared hash is ~50% noise on a novel species, 2 hashes rarely. Containment above ~0.1% costs real
# bacterial contigs (large contigs share <1% of hashes with the nearest GTDB genome), but 0.1% stops
# chromosome-scale scaffolds being called on a few dozen hashes. Margin 2 is flat-to-better on F1
# than 1 and more conservative.
DEFAULT_MIN_HASHES = 2.0
DEFAULT_MIN_CONTAINMENT = 0.001
DEFAULT_MARGIN = 2.0


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
    """Return the trailing integer of a contig name (e.g. 'scaffold_42' -> 42), or inf if none."""
    match = re.search(r'(\d+)$', query_name)
    return int(match.group(1)) if match else float('inf')


# ---------------------------------------------------------------------------
# Exclusions
# ---------------------------------------------------------------------------

def read_exclude_file(path):
    """Read accessions to exclude from a file.

    Accepts a TSV/CSV with an 'assembly_accession' column (other columns ignored), or a plain list
    with one accession per line (first whitespace/comma/tab-separated token). Blank lines and lines
    starting with '#' are skipped.
    """
    validate_file_paths([path], "Exclude accessions")
    with open(path) as fh:
        lines = [ln.rstrip('\n') for ln in fh if ln.strip() and not ln.lstrip().startswith('#')]
    if not lines:
        return set()
    first = lines[0]
    sep = '\t' if '\t' in first else (',' if ',' in first else None)
    if sep is not None:
        header = [c.strip().strip('"') for c in first.split(sep)]
        if 'assembly_accession' in header:
            idx = header.index('assembly_accession')
            accs = set()
            for ln in lines[1:]:
                cols = ln.split(sep)
                if idx < len(cols) and cols[idx].strip().strip('"'):
                    accs.add(cols[idx].strip().strip('"'))
            return accs
    return {re.split(r'[\s,]+', ln.strip())[0] for ln in lines}


def tolid_accessions(assembly_df, tolid):
    """Accessions in the taxa table whose assembly_name belongs to this ToLID (same specimen).

    assembly_name looks like 'icCatNigr1.hap1.1', 'icCatNigr1.1', 'icCatNigr1_alternate_haplotype'.
    Matches '<tolid>' followed by '.', '_' or end of string, so icCatNigr1 does not match icCatNigr10.
    """
    if not tolid:
        return set()
    if not TOLID_RE.match(tolid):
        logger.warning(f"'{tolid}' does not look like a ToLID; own-specimen exclusion skipped")
        return set()
    if assembly_df is None or 'assembly_name' not in assembly_df.columns:
        logger.warning("Assembly database has no 'assembly_name' column; own-specimen exclusion skipped")
        return set()
    pattern = f"^{re.escape(tolid)}([._]|$)"
    accs = set(
        assembly_df
        .filter(pl.col('assembly_name').cast(pl.Utf8).str.contains(pattern))
        .get_column('assembly_accession')
        .to_list()
    )
    logger.info(f"Own-specimen exclusion ({tolid}): {len(accs)} accession(s) {', '.join(sorted(accs))}")
    return accs


def build_exclude_set(exclude_accessions=None, exclude_files=None, assembly_df=None, exclude_tolid=None):
    """Combine CLI accessions (comma-separated allowed), exclude files and own-ToLID accessions."""
    exclude_set = set()
    for acc_arg in exclude_accessions or []:
        for acc in acc_arg.split(','):
            if acc.strip():
                exclude_set.add(acc.strip())
    for path in exclude_files or []:
        accs = read_exclude_file(path)
        logger.info(f"Read {len(accs)} accessions to exclude from {path}")
        exclude_set |= accs
    exclude_set |= tolid_accessions(assembly_df, exclude_tolid)
    if exclude_set:
        shown = ', '.join(sorted(exclude_set)[:20]) + (' ...' if len(exclude_set) > 20 else '')
        logger.info(f"Excluding {len(exclude_set)} accessions: {shown}")
    return exclude_set


def exclusion_mask(accession_col, exclude_set):
    """Polars expression: True where accession is excluded.

    Entries with a version ('GCA_1.2') match exactly; entries without one ('GCA_1') match any version.
    """
    versioned = [a for a in exclude_set if '.' in a]
    unversioned = [a for a in exclude_set if '.' not in a]
    col = pl.col(accession_col)
    expr = col.is_in(versioned)
    if unversioned:
        expr = expr | col.str.replace(r'\.\d+$', '').is_in(unversioned)
    return expr


# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------

def db_label(file_path):
    """DB name from a multisearch file name: '<sample>_vs_<db>.multisearch_results.csv' -> '<db>'."""
    name = os.path.basename(file_path)
    name = re.sub(r'^.*_vs_', '', name)
    return re.sub(r'\.multisearch_results(\.csv)?$', '', name)


def parse_sourmash_results(file_paths, exclude_accessions=None, exclude_set=None):
    """Read multisearch CSVs into one polars DataFrame of hits.

    Columns: query_name, match_name (accession only), containment, jaccard, intersect_hashes,
    ksize, scaled, source. Duplicate header rows and non-numeric rows are dropped.
    Missing ksize/scaled columns (old/hand-made inputs) default to null / 200.

    Returns (hits_df, all_queries) where all_queries is every contig name seen before exclusions,
    so contigs whose only hits were excluded can still be reported as 'no_hit'.
    """
    validate_file_paths(file_paths, "Sourmash results")
    logger.info(f"Processing {len(file_paths)} sourmash results file(s)")

    all_dfs = []
    for file_path in file_paths:
        try:
            df = pl.read_csv(file_path, infer_schema_length=0)  # all Utf8: duplicate headers are possible
        except pl.exceptions.NoDataError:
            logger.warning(f"Sourmash results file {file_path} is empty or has no columns")
            continue
        if df.is_empty():
            logger.warning(f"Sourmash results file {file_path} is empty")
            continue
        missing = [c for c in ['query_name', 'match_name', 'containment', 'jaccard'] if c not in df.columns]
        if missing:
            raise ValueError(f"Missing required columns in {file_path}: {missing}")
        keep = [c for c in ['query_name', 'match_name', 'containment', 'jaccard', 'intersect_hashes', 'ksize', 'scaled'] if c in df.columns]
        df = df.select(keep).with_columns(pl.lit(db_label(file_path)).alias('source'),
                                          pl.lit('_vs_' in os.path.basename(file_path)).alias('source_is_db'))
        logger.info(f"Loaded {len(df)} sourmash results from {file_path}")
        all_dfs.append(df)

    if not all_dfs:
        raise ValueError("No valid sourmash results files provided")

    hits = pl.concat(all_dfs, how="diagonal_relaxed")
    for c in ['intersect_hashes', 'ksize', 'scaled']:
        if c not in hits.columns:
            hits = hits.with_columns(pl.lit(None, dtype=pl.Utf8).alias(c))

    hits = (
        hits
        .filter(pl.col("containment").str.contains(r"^[0-9]") & pl.col("jaccard").str.contains(r"^[0-9]"))
        .with_columns(
            pl.col("containment").cast(pl.Float64),
            pl.col("jaccard").cast(pl.Float64),
            pl.col("intersect_hashes").cast(pl.Float64, strict=False).fill_null(0.0),
            pl.col("ksize").cast(pl.Float64, strict=False).cast(pl.Int64, strict=False),
            pl.col("scaled").cast(pl.Float64, strict=False).fill_null(float(REFERENCE_SCALED)).cast(pl.Int64),
            # match_name is '<accession> <description>' in multisearch output
            pl.col("match_name").str.split(" ").list.first(),
            # query_name is the full FASTA header; keep the sequence id (first token), as FCS-GX and the
            # .genome / .fai contig lists do
            pl.col("query_name").str.strip_chars().str.split(" ").list.first(),
        )
    )
    logger.info(f"Combined sourmash results: {len(hits)} valid rows")

    all_queries = hits.get_column('query_name').unique().to_list()

    if exclude_set is None:
        exclude_set = build_exclude_set(exclude_accessions)
    if exclude_set:
        before = len(hits)
        hits = hits.filter(~exclusion_mask('match_name', exclude_set))
        logger.info(f"Filtered out {before - len(hits)} matches due to excluded accessions")
    logger.info(f"Total matches after filtering: {len(hits)}")
    return hits, all_queries


def parse_assembly_database(file_paths):
    """Parse the assembly (taxa) database files into one polars DataFrame, deduplicated on accession."""
    validate_file_paths(file_paths, "Assembly database")
    logger.info(f"Processing {len(file_paths)} assembly database file(s)")

    all_dfs = []
    for file_path in file_paths:
        try:
            assembly_df = pl.read_csv(file_path, infer_schema_length=0)  # all Utf8: lineage columns are text
        except pl.exceptions.NoDataError:
            logger.warning(f"Assembly database file {file_path} is empty or has no columns")
            continue
        if assembly_df.is_empty():
            logger.warning(f"Assembly database file {file_path} is empty")
            continue
        logger.info(f"Loaded {len(assembly_df)} assembly records from {file_path}")
        all_dfs.append(assembly_df)

    if not all_dfs:
        raise ValueError("No valid assembly database files provided")

    combined_df = pl.concat(all_dfs, how="diagonal_relaxed")
    logger.info(f"Combined assembly database: {len(combined_df)} total rows")

    missing_required = [c for c in ['assembly_accession', 'taxid'] if c not in combined_df.columns]
    if missing_required:
        raise ValueError(
            f"Assembly database is missing required column(s): {missing_required}. "
            f"Found columns: {combined_df.columns}"
        )

    before = len(combined_df)
    combined_df = combined_df.unique(subset=['assembly_accession'], keep="first", maintain_order=True)
    if len(combined_df) < before:
        logger.warning("Duplicate 'assembly_accession' values found after merging. Keeping first occurrence.")

    initial_count = len(combined_df)
    combined_df = combined_df.with_columns(
        pl.col('taxid').str.strip_chars().str.replace(r'\.0$', '')
    ).filter(
        pl.col('taxid').is_not_null() & ~pl.col('taxid').is_in(['nan', ''])
    )
    logger.info(f"Filtered out {initial_count - len(combined_df)} assemblies without valid taxid")
    return combined_df


def get_target_genomes(assembly_df, target_taxa):
    """Set of accessions matching the target taxa.

    Returns:
      - empty set  → target_taxa not provided (no target/non-target classification)
      - None       → target_taxa provided but taxon absent from DB (conservative: everything is target)
      - non-empty  → normal case
    """
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
        level_matches = (
            assembly_df
            .filter(pl.col(level).cast(pl.Utf8).str.to_lowercase().str.strip_chars() == target_value)
            .get_column('assembly_accession')
            .to_list()
        )
        target_genomes.update(level_matches)
        logger.info(f"Found {len(level_matches)} genomes matching {level}={target_value}")

    if not target_genomes:
        # TODO: automatic taxonomic level fallback (order → class → phylum) for rare taxa.
        logger.warning(
            f"Target taxa {target_taxa} not found in assembly database — no genomes matched. "
            f"Applying conservative fallback: all hits are treated as is_target=True to prevent "
            f"accidental removal. Consider a higher taxonomic level (class/phylum)."
        )
        return None

    logger.info(f"Total target genomes: {len(target_genomes)}")
    return target_genomes


# ---------------------------------------------------------------------------
# Scoring and classification
# ---------------------------------------------------------------------------

def resolve_target_taxon(assembly_df, target_taxa, exclude_set=None, min_genomes=1):
    """Pick the target taxon from an ordered preference list (e.g. order, class, phylum, kingdom).

    Returns ({level: value}, n_genomes) for the first level whose taxon has at least `min_genomes`
    DB genomes left after the exclusions (own ToLID / species / exclude file), or (None, 0) if none does.
    """
    exclude_set = exclude_set or set()
    for level, value in (target_taxa or {}).items():
        if assembly_df is None or level not in assembly_df.columns:
            logger.warning(f"Target level '{level}' not in the assembly database columns; skipped")
            continue
        accs = (assembly_df
                .filter(pl.col(level).cast(pl.Utf8).str.to_lowercase().str.strip_chars() == value)
                .filter(~exclusion_mask('assembly_accession', exclude_set) if exclude_set else pl.lit(True))
                .height)
        if accs >= min_genomes:
            logger.info(f"Target taxon: {level}:{value} ({accs} DB genomes after exclusions, minimum {min_genomes})")
            return {level: value}, accs
        logger.warning(f"Target {level}:{value} has {accs} DB genomes after exclusions (< {min_genomes}); "
                       f"falling back to the next level")
    return None, 0


def annotate_hits(hits, assembly_df=None, target_taxa=None):
    """Add taxa, db, is_target, norm_hashes, score (= containment) and per-contig rank (top_n).

    is_target: True / False, or null when no target taxa were given. Accessions missing from the
    taxa table are treated as target (never remove on unknown evidence).
    """
    target_genomes = get_target_genomes(assembly_df, target_taxa)

    if assembly_df is not None:
        extra = [c for c in ['sourmash_db_lineage'] if c in assembly_df.columns]
        lookup = assembly_df.select(
            [pl.col('assembly_accession').alias('match_name'), pl.col('taxid').alias('taxa')]
            + [pl.col(c) for c in extra]
        )
        hits = hits.join(lookup, on='match_name', how='left')
        if 'sourmash_db_lineage' not in hits.columns:
            hits = hits.with_columns(pl.lit(None, dtype=pl.Utf8).alias('sourmash_db_lineage'))
    else:
        hits = hits.with_columns(pl.lit(None, dtype=pl.Utf8).alias('taxa'),
                                 pl.lit(None, dtype=pl.Utf8).alias('sourmash_db_lineage'))

    unknown = pl.col('taxa').is_null()
    n_unknown = hits.filter(unknown).get_column('match_name').n_unique()
    if n_unknown and assembly_df is not None:
        logger.warning(f"{n_unknown} matched accession(s) not found in assembly_taxa database — treated as target (safe side)")

    if target_genomes is None:
        is_target = pl.lit(True)
    elif target_genomes:
        is_target = unknown | pl.col('match_name').is_in(list(target_genomes))
    else:
        is_target = pl.lit(None, dtype=pl.Boolean)

    hits = hits.with_columns(
        pl.col('taxa').fill_null('nan'),
        # file name when it follows the SOURMASH_MULTISEARCH '<sample>_vs_<db>' naming, else the taxa table's DB
        pl.when(pl.col('source_is_db')).then(pl.col('source'))
        .otherwise(pl.coalesce(pl.col('sourmash_db_lineage'), pl.col('source'))).alias('db'),
        is_target.alias('is_target'),
        (pl.col('intersect_hashes') * pl.col('scaled') / REFERENCE_SCALED).alias('norm_hashes'),
        # ranking score. Containment is scale-invariant; containment-ANI (C^(1/k)) was tried as a
        # k-invariant alternative and did no better on the F1 benchmark (TODO 04)
        pl.col('containment').alias('score'),
    ).drop('sourmash_db_lineage', 'source_is_db')

    hits = hits.sort(['query_name', 'score', 'norm_hashes', 'match_name'],
                     descending=[False, True, True, False])
    hits = hits.with_columns(
        pl.int_range(1, pl.len() + 1).over('query_name').alias('top_n')
    )
    return hits


def evidence_expr(min_hashes, min_containment):
    """Polars expression: hit carries enough evidence to support a non-target call."""
    return (pl.col('norm_hashes') >= min_hashes) & (pl.col('containment') >= min_containment)


def classify_contigs(hits, all_queries, min_hashes=DEFAULT_MIN_HASHES,
                     min_containment=DEFAULT_MIN_CONTAINMENT, margin=DEFAULT_MARGIN):
    """One row per contig: FCS-GX-style action plus best target / best evidence non-target / top hit.

    action (FCS-GX vocabulary):
      EXCLUDE: best evidence non-target containment > margin * best target containment (or no target hit)
      REVIEW : best evidence non-target beats the best target hit, but not by the margin (ambiguous)
      NA     : everything else — target wins, only weak non-target evidence, no hits, or no target taxa.
    """
    ev = evidence_expr(min_hashes, min_containment)
    per_contig = (
        hits
        .with_columns(ev.alias('_ev'))
        .group_by('query_name')
        .agg(
            pl.len().alias('n_hits'),
            pl.col('match_name').first().alias('top_accession'),
            pl.col('score').filter(pl.col('is_target') == True).max().alias('best_target_score'),  # noqa: E712
            pl.col('match_name').filter(pl.col('is_target') == True).first().alias('best_target_accession'),  # noqa: E712
            pl.col('score').filter((pl.col('is_target') == False) & pl.col('_ev')).max().alias('best_non_target_score'),  # noqa: E712
            pl.col('match_name').filter((pl.col('is_target') == False) & pl.col('_ev')).first().alias('best_non_target_accession'),  # noqa: E712
            pl.col('is_target').is_null().all().alias('_no_taxa'),
        )
    )
    # hits are sorted by score desc within contig, so .first() above is the best hit
    queries = pl.DataFrame({'query_name': sorted(set(all_queries), key=lambda q: (extract_contig_number(q), q))},
                           schema={'query_name': pl.Utf8})
    df = queries.join(per_contig, on='query_name', how='left', maintain_order='left')

    t = pl.col('best_target_score')
    nt = pl.col('best_non_target_score')
    exclude = t.is_null() | (nt > (pl.lit(margin, dtype=pl.Float64) * t)) if not math.isinf(margin) else t.is_null()
    df = df.with_columns(
        pl.when(pl.col('n_hits').is_null() | pl.col('_no_taxa') | nt.is_null()).then(pl.lit('NA'))
        .when(exclude).then(pl.lit('EXCLUDE'))
        .when(nt > t).then(pl.lit('REVIEW'))
        .otherwise(pl.lit('NA'))
        .alias('action'),
        pl.col('n_hits').fill_null(0),
    ).drop('_no_taxa')
    # the hit reported for the contig: the non-target evidence for EXCLUDE/REVIEW, else the best target
    # hit, else the top hit
    df = df.with_columns(
        pl.when(pl.col('action') != 'NA').then(pl.col('best_non_target_accession'))
        .otherwise(pl.coalesce(pl.col('best_target_accession'), pl.col('top_accession')))
        .alias('reported_accession')
    )
    return df


# ---------------------------------------------------------------------------
# NCBI taxonomy: FCS-GX divisions (gx-div) from BLAST names, species names
# ---------------------------------------------------------------------------

def _read_dmp(path, columns):
    """Read an NCBI *.dmp file ('\t|\t'-separated) keeping the given 0-based field indices."""
    df = pl.read_csv(path, separator='\t', has_header=False, quote_char=None, infer_schema_length=0,
                     truncate_ragged_lines=True)
    # fields sit at even column positions (odd ones are the '|' separators)
    return df.select([pl.col(df.columns[2 * i]).alias(name) for i, name in columns.items()])


def load_gx_divs(taxdump_dir, mapping_file, taxids):
    """taxid -> FCS-GX division (e.g. 'prok:CFB group bacteria', 'anml:insects').

    FCS-GX derives its gx-div from the NCBI BLAST name of the taxon (nearest ancestor with a
    'blast name' in names.dmp) mapped through the BLAST-div → GX-div table (make_gxdb
    blast_names_mapping.tsv). Taxids without a BLAST name in their lineage get 'unkn:unknown'.
    """
    taxids = {str(t) for t in taxids if t is not None and str(t).isdigit()}
    if not taxids:
        return {}
    missing = [f for f in ('names.dmp', 'taxidlineage.dmp') if not os.path.exists(os.path.join(taxdump_dir, f))]
    if missing:
        logger.warning(f"{taxdump_dir} has no {', '.join(missing)}; sourmash_div left empty")
        return {}
    mapping = {}
    with open(mapping_file) as fh:
        for ln in fh:
            if ln.startswith('#') or not ln.strip():
                continue
            blast_div, gx_div = ln.rstrip('\n').split('\t')[:2]
            mapping[blast_div] = gx_div
    names = _read_dmp(os.path.join(taxdump_dir, 'names.dmp'), {0: 'taxid', 1: 'name', 3: 'class'})
    blast_name = dict(names.filter(pl.col('class') == 'blast name').select('taxid', 'name').iter_rows())
    merged = {}
    merged_path = os.path.join(taxdump_dir, 'merged.dmp')
    if os.path.exists(merged_path) and os.path.getsize(merged_path) > 0:
        merged = dict(_read_dmp(merged_path, {0: 'old', 1: 'new'}).filter(pl.col('old').is_in(list(taxids))).iter_rows())
    current = {t: merged.get(t, t) for t in taxids}
    lineages = dict(
        _read_dmp(os.path.join(taxdump_dir, 'taxidlineage.dmp'), {0: 'taxid', 1: 'lineage'})
        .filter(pl.col('taxid').is_in(list(set(current.values()))))
        .iter_rows()
    )
    divs = {}
    for t, cur in current.items():
        chain = [cur] + list(reversed((lineages.get(cur) or '').split()))
        bn = next((blast_name[x] for x in chain if x in blast_name), None)
        divs[t] = mapping.get(bn, 'unkn:unknown') if bn else 'unkn:unknown'
    logger.info(f"gx-div assigned for {len(divs)} taxids ({sum(v == 'unkn:unknown' for v in divs.values())} unknown)")
    return divs


def species_name_for_taxid(taxid, taxdump_dir=None, assembly_df=None):
    """Species name of a taxid: rankedlineage.dmp (species column, or name if the taxid is a species),
    else the taxa table rows with that taxid."""
    taxid = str(taxid)
    if taxdump_dir and os.path.exists(os.path.join(taxdump_dir, 'rankedlineage.dmp')):
        rl = (_read_dmp(os.path.join(taxdump_dir, 'rankedlineage.dmp'), {0: 'taxid', 1: 'name', 2: 'species'})
              .filter(pl.col('taxid') == taxid))
        if rl.height:
            name, species = rl.row(0)[1], rl.row(0)[2]
            return (species or name or '').strip() or None
    if assembly_df is not None:
        rows = assembly_df.filter(pl.col('taxid') == taxid)
        for col in ['species', 'tax_name']:
            if col in rows.columns:
                vals = [v for v in rows.get_column(col).to_list() if v and v.strip() not in ('', 'nan')]
                if vals:
                    return vals[0].strip()
    return None


def species_accessions(assembly_df, species, taxdump_dir=None):
    """Accessions of DB assemblies of the given species (a name, or a taxid resolved to its species).

    Matches taxid (when a taxid is given), or the 'species' or 'tax_name' column equal to the species
    name (case-insensitive). Subspecies/strain rows match through their 'species' column."""
    if not species or assembly_df is None:
        return set()
    taxid = str(species).strip() if str(species).strip().isdigit() else None
    name = species_name_for_taxid(taxid, taxdump_dir, assembly_df) if taxid else str(species).strip()
    cond = pl.lit(False)
    if taxid:
        cond = cond | (pl.col('taxid') == taxid)
    if name:
        low = name.lower()
        for col in ['species', 'tax_name']:
            if col in assembly_df.columns:
                # exact match only: tax_name prefixes also catch e.g. 'Bathymodiolus heckerae thiotrophic
                # gill symbiont' (a contaminant-like symbiont genome named after its host)
                c = pl.col(col).cast(pl.Utf8).str.to_lowercase().str.strip_chars()
                cond = cond | (c == low)
    accs = set(assembly_df.filter(cond).get_column('assembly_accession').to_list())
    logger.info(f"Own-species exclusion ({species} → {name}): {len(accs)} accession(s) {', '.join(sorted(accs)[:20])}")
    return accs


def read_contig_sizes(path):
    """Contig names and lengths from a sizes / .genome / .fai file (first two tab-separated columns)."""
    validate_file_paths([path], "Contig sizes")
    df = pl.read_csv(path, separator='\t', has_header=False, infer_schema_length=0, truncate_ragged_lines=True)
    return df.select(pl.col(df.columns[0]).alias('query_name'),
                     pl.col(df.columns[1]).cast(pl.Int64, strict=False).alias('seq_len'))


def _is_target_str():
    """is_target as 'True' / 'False' / 'None' strings (parser 1.x format)."""
    return (pl.when(pl.col('is_target').is_null()).then(pl.lit('None'))
            .when(pl.col('is_target')).then(pl.lit('True'))
            .otherwise(pl.lit('False')).alias('is_target'))


def _contig_order(df, col='query_name'):
    """Sort key per contig: trailing contig number, then name."""
    names = df.get_column(col).unique().to_list()
    return pl.DataFrame({col: names, '_n': [float(extract_contig_number(n)) for n in names]},
                        schema={col: pl.Utf8, '_n': pl.Float64})


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_summary_output(hits, output_file, min_hashes=DEFAULT_MIN_HASHES, min_containment=DEFAULT_MIN_CONTAINMENT):
    """One row per hit. The first 8 columns are unchanged from parser 1.x (top_n now ranks by score)."""
    order = _contig_order(hits)
    out = (
        hits
        .with_columns(evidence_expr(min_hashes, min_containment).alias('passes_thresholds'))
        .join(order, on='query_name', how='left')
        .sort(['_n', 'query_name', 'top_n'])
        .select(
            pl.col('query_name').alias('header'),
            pl.col('match_name').alias('assembly_accession'),
            'taxa',
            'top_n',
            pl.col('containment').round(6),
            pl.col('jaccard').round(9),
            'intersect_hashes',
            _is_target_str(),
            'db', 'ksize', 'scaled',
            pl.col('norm_hashes').round(3),
            'passes_thresholds',
        )
    )
    out.write_csv(output_file)
    logger.info(f"Wrote {len(out)} hit rows to {output_file}")


def lineage_table(assembly_df):
    """accession -> lineage columns; empty species falls back to tax_name (e.g. 'uncultured virus')."""
    cols = [c for c in LINEAGE_COLUMNS if c in assembly_df.columns]
    df = assembly_df.select(['assembly_accession'] + cols + (['tax_name'] if 'tax_name' in assembly_df.columns else []))
    for c in LINEAGE_COLUMNS:
        if c not in df.columns:
            df = df.with_columns(pl.lit(None, dtype=pl.Utf8).alias(c))
        df = df.with_columns(
            pl.when(pl.col(c).is_null() | pl.col(c).str.strip_chars().is_in(['', 'nan', 'None']))
            .then(None).otherwise(pl.col(c).str.strip_chars()).alias(c)
        )
    if 'tax_name' in df.columns:
        df = df.with_columns(pl.coalesce(pl.col('species'), pl.col('tax_name')).alias('species')).drop('tax_name')
    return df.rename({'assembly_accession': 'match_name'})


def write_non_target_output(classified, hits, assembly_df, output_file):
    """Contigs with action EXCLUDE, with their best evidence non-target hit and its lineage.

    Kept for backward compatibility (generate_contamination_bed.py, older autofilter inputs).
    Columns 1-13 are unchanged from parser 1.x; extra columns are appended.
    """
    nt = classified.filter(pl.col('action') == 'EXCLUDE').select(
        'query_name', pl.col('best_non_target_accession').alias('match_name'), 'best_target_score', 'best_target_accession')
    rows = nt.join(hits, on=['query_name', 'match_name'], how='left', maintain_order='left')
    rows = rows.unique(subset=['query_name'], keep='first', maintain_order=True)
    rows = rows.join(lineage_table(assembly_df), on='match_name', how='left', maintain_order='left')
    out = rows.select(
        pl.col('query_name').alias('header'),
        pl.col('match_name').alias('assembly_accession'),
        'taxa',
        pl.col('containment').round(6),
        pl.col('jaccard').round(9),
        'intersect_hashes',
        *[pl.col(c).fill_null('None') for c in LINEAGE_COLUMNS],
        'db', 'ksize', 'scaled',
        pl.col('norm_hashes').round(3),
        pl.col('best_target_accession'),
        pl.col('best_target_score').round(6).alias('best_target_containment'),
    )
    out.write_csv(output_file)
    logger.info(f"Non-target (EXCLUDE) contigs: {len(out)} → {output_file}")


# Column layout mirrors ascc's parsed FCS-GX table (parse_fcsgx_result.py: scaff, fcs_gx_top_tax_name,
# fcs_gx_top_taxid, fcs_gx_div, ..., fcs_gx_action, ..., fcs_gx_species ... fcs_gx_domain), then
# sourmash-specific evidence for the reported hit.
CLASSIFICATION_COLUMNS = [
    'scaff', 'sourmash_top_tax_name', 'sourmash_top_taxid', 'sourmash_div', 'sourmash_action',
    'sourmash_species', 'sourmash_genus', 'sourmash_family', 'sourmash_order', 'sourmash_class',
    'sourmash_phylum', 'sourmash_kingdom', 'sourmash_domain',
    'sourmash_seq_len', 'sourmash_top_accession', 'sourmash_top_is_target', 'sourmash_db', 'sourmash_ksize',
    'sourmash_scaled', 'sourmash_containment', 'sourmash_intersect_hashes', 'sourmash_norm_hashes',
    'sourmash_best_target_accession', 'sourmash_best_target_containment', 'sourmash_n_hits',
    'sourmash_target_taxon',
]


def build_classification_table(classified, hits, assembly_df=None, gx_divs=None, sizes=None, target_used=None):
    """One row per contig, FCS-GX-style (see CLASSIFICATION_COLUMNS)."""
    rep = classified.select('query_name', 'action', 'n_hits', 'best_target_accession', 'best_target_score',
                            pl.col('reported_accession').alias('match_name'))
    rep = rep.join(hits.select('query_name', 'match_name', 'taxa', 'is_target', 'db', 'ksize', 'scaled',
                               'containment', 'intersect_hashes', 'norm_hashes'),
                   on=['query_name', 'match_name'], how='left', maintain_order='left')
    rep = rep.unique(subset=['query_name'], keep='first', maintain_order=True)
    if assembly_df is not None:
        lin = lineage_table(assembly_df)
        extra = [c for c in ['tax_name', 'superkingdom'] if c in assembly_df.columns]
        if extra:
            lin = lin.join(assembly_df.select(pl.col('assembly_accession').alias('match_name'), *extra), on='match_name', how='left')
        rep = rep.join(lin, on='match_name', how='left', maintain_order='left')
    for c in LINEAGE_COLUMNS + ['tax_name', 'superkingdom']:
        if c not in rep.columns:
            rep = rep.with_columns(pl.lit(None, dtype=pl.Utf8).alias(c))
    divs = gx_divs or {}
    rep = rep.with_columns(
        pl.col('taxa').replace('nan', None),
        pl.col('taxa').replace_strict(divs, default=None, return_dtype=pl.Utf8).alias('div'),
    )
    if sizes is not None:
        rep = rep.join(sizes, on='query_name', how='left', maintain_order='left')
    else:
        rep = rep.with_columns(pl.lit(None, dtype=pl.Int64).alias('seq_len'))
    out = rep.select(
        pl.col('query_name').alias('scaff'),
        pl.coalesce(pl.col('tax_name'), pl.col('species')).alias('sourmash_top_tax_name'),
        pl.col('taxa').alias('sourmash_top_taxid'),
        pl.col('div').alias('sourmash_div'),
        pl.col('action').alias('sourmash_action'),
        *[pl.col(c).alias(f'sourmash_{c}') for c in LINEAGE_COLUMNS],
        pl.col('superkingdom').alias('sourmash_domain'),
        pl.col('seq_len').alias('sourmash_seq_len'),
        pl.col('match_name').alias('sourmash_top_accession'),
        pl.col('is_target').alias('sourmash_top_is_target'),
        pl.col('db').alias('sourmash_db'),
        pl.col('ksize').alias('sourmash_ksize'),
        pl.col('scaled').alias('sourmash_scaled'),
        pl.col('containment').round(6).alias('sourmash_containment'),
        pl.col('intersect_hashes').alias('sourmash_intersect_hashes'),
        pl.col('norm_hashes').round(3).alias('sourmash_norm_hashes'),
        pl.col('best_target_accession').alias('sourmash_best_target_accession'),
        pl.col('best_target_score').round(6).alias('sourmash_best_target_containment'),
        pl.col('n_hits').alias('sourmash_n_hits'),
        # target taxon actually used (after any fallback); 'none' = no level in the DB, all hits target
        pl.lit(target_used, dtype=pl.Utf8).alias('sourmash_target_taxon'),
    )
    return out


def write_classification_output(table, output_file):
    table.write_csv(output_file)
    counts = dict(table.group_by('sourmash_action').len().iter_rows())
    logger.info(f"Actions: {counts} → {output_file}")


def output_basename(sourmash_files, target_taxa, prefix=None, first_only=False):
    if prefix:
        base = prefix
    else:
        base = "sourmash_results" if len(sourmash_files) > 1 else os.path.splitext(os.path.basename(sourmash_files[0]))[0]
    items = list((target_taxa or {}).items())[:1] if first_only else (target_taxa or {}).items()
    taxa_suffix = "".join(f"_{level}_{value}" for level, value in items)
    return f"{base}{taxa_suffix}"


def main(sourmash_files=None, assembly_dbs=None, target_taxa=None, outdir=None, exclude_accessions=None,
         exclude_files=None, exclude_tolid=None, exclude_species=None, min_hashes=DEFAULT_MIN_HASHES,
         min_containment=DEFAULT_MIN_CONTAINMENT, margin=DEFAULT_MARGIN, prefix=None,
         contig_sizes=None, taxdump=None, div_mapping=None, target_fallback=False, min_target_genomes=1):
    """Run the parser; returns the per-contig classification table (FCS-GX-style columns).

    target_fallback: target_taxa is an ordered preference list (e.g. order, class, phylum, kingdom); the
    first taxon with >= min_target_genomes DB genomes (after exclusions) is used. Without it, several
    target taxa are a union (parser <= 2.1 behaviour).
    """
    logger.info("=" * 60)
    logger.info(f"sourmash_taxonomy_parser {__version__}")
    logger.info(f"  Sourmash results files: {', '.join(sourmash_files)}")
    logger.info(f"  Assembly database files: {', '.join(assembly_dbs) if assembly_dbs else 'N/A'}")
    logger.info(f"  Target taxa: {target_taxa or 'None'}")
    logger.info(f"  min_hashes(s{REFERENCE_SCALED})={min_hashes} min_containment={min_containment} margin={margin}")
    logger.info("=" * 60)

    assembly_df = parse_assembly_database(assembly_dbs) if assembly_dbs else None
    exclude_set = build_exclude_set(exclude_accessions, exclude_files, assembly_df, exclude_tolid)
    exclude_set |= species_accessions(assembly_df, exclude_species, taxdump)
    hits, all_queries = parse_sourmash_results(sourmash_files, exclude_set=exclude_set)
    sizes = None
    if contig_sizes:
        sizes = read_contig_sizes(contig_sizes)
        only_hits = set(all_queries) - set(sizes.get_column('query_name').to_list())
        if only_hits:
            logger.warning(f"{len(only_hits)} contig(s) in the multisearch results are not in {contig_sizes}")
        all_queries = list(set(all_queries) | set(sizes.get_column('query_name').to_list()))
    requested = target_taxa
    if target_fallback and target_taxa:
        chosen, _ = resolve_target_taxon(assembly_df, target_taxa, exclude_set, min_target_genomes)
        # no level qualifies: keep the requested (first) taxon -> conservative 'everything is target'
        # if it has no genomes at all
        target_taxa = chosen or dict(list(target_taxa.items())[:1])
        target_used = (':'.join(next(iter(chosen.items()))) if chosen else 'none')
    else:
        target_used = ' '.join(f"{k}:{v}" for k, v in (target_taxa or {}).items()) or None
    if target_taxa and get_target_genomes(assembly_df, target_taxa) is None:
        target_used = 'none'
    logger.info(f"Target taxon used: {target_used}")
    hits = annotate_hits(hits, assembly_df, target_taxa)
    classified = classify_contigs(hits, all_queries, min_hashes, min_containment, margin)

    gx_divs = {}
    if taxdump and div_mapping:
        needed = hits.filter(pl.col('match_name').is_in(classified.get_column('reported_accession').drop_nulls().to_list()))
        gx_divs = load_gx_divs(taxdump, div_mapping, needed.get_column('taxa').unique().to_list())
    elif taxdump or div_mapping:
        logger.warning("sourmash_div needs both --taxdump and --div_mapping; leaving it empty")
    table = build_classification_table(classified, hits, assembly_df, gx_divs, sizes, target_used)

    os.makedirs(outdir, exist_ok=True)
    # file names keep the requested taxon (first level), whatever level was used
    base = os.path.join(outdir, output_basename(sourmash_files, requested, prefix, first_only=target_fallback))

    write_summary_output(hits, f"{base}.summary.csv", min_hashes, min_containment)
    write_classification_output(table, f"{base}.classification.csv")
    if assembly_df is not None and target_taxa:
        write_non_target_output(classified, hits, assembly_df, f"{base}.non_target.csv")
    return table


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Parse sourmash multisearch results into per-hit / per-contig target vs non-target calls.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--sourmash_results', nargs='+', help='Path(s) to the sourmash multisearch results file(s)')
    parser.add_argument('-a', '--assembly_db', nargs='+', help='Path(s) to assembly database file(s) with taxonomic information')
    parser.add_argument('--target_taxa', nargs='+', help='Target taxa in format taxon:value (e.g., order:coleoptera)')
    parser.add_argument('--target_fallback', action='store_true',
                        help='Treat --target_taxa as an ordered list (e.g. order:x class:y phylum:z kingdom:w) and use the first '
                             'taxon with at least --min_target_genomes DB genomes after exclusions')
    parser.add_argument('--min_target_genomes', type=int, default=1,
                        help='With --target_fallback: minimum DB genomes of the target taxon')
    parser.add_argument('--accessions_to_exclude', nargs='+', default=None,
                        help='Assembly accessions to drop from the hits (comma-separated or multiple values). '
                             'An accession without version matches every version.')
    parser.add_argument('--exclude_accessions_file', nargs='+', default=None,
                        help='File(s) of accessions to drop: TSV/CSV with an assembly_accession column, or one accession per line')
    parser.add_argument('--exclude_tolid', default=None,
                        help="Drop every DB assembly of this specimen: taxa-table assembly_name starting with '<tolid>.' or '<tolid>_'")
    parser.add_argument('--exclude_species', default=None,
                        help="Drop every DB assembly of this species: a taxid (resolved with --taxdump rankedlineage.dmp, "
                             "else the taxa table) or a species name. For testing on published assemblies only.")
    parser.add_argument('--contig_sizes', default=None,
                        help='Sizes/.genome/.fai file (name<TAB>length): every contig gets a row in classification.csv')
    parser.add_argument('--taxdump', default=None,
                        help='NCBI new_taxdump directory (names.dmp, taxidlineage.dmp, merged.dmp, rankedlineage.dmp) '
                             'for sourmash_div and --exclude_species taxids')
    parser.add_argument('--div_mapping', default=None,
                        help='BLAST-div → FCS-GX div table (e.g. assets/fcs_gx_blast_div_mapping.tsv)')
    parser.add_argument('--min_hashes', type=float, default=DEFAULT_MIN_HASHES,
                        help=f'Minimum shared hashes, normalised to scaled={REFERENCE_SCALED} (intersect_hashes*scaled/{REFERENCE_SCALED}), '
                             'for a non-target hit to count as evidence')
    parser.add_argument('--min_containment', type=float, default=DEFAULT_MIN_CONTAINMENT,
                        help='Minimum containment (fraction of contig hashes) for a non-target hit to count as evidence')
    parser.add_argument('--margin', type=float, default=DEFAULT_MARGIN,
                        help="Call non-target only if best evidence non-target containment > margin * best target containment. "
                             "'inf' = any target hit keeps the contig (parser 1.x behaviour)")
    parser.add_argument('--log', help='Path to log file (if not provided, logs to stderr)', default=None)
    parser.add_argument('-o', '--outdir', help='Output directory for results', default=None)
    parser.add_argument('--prefix', default=None,
                        help="Output file prefix (default: input basename, or 'sourmash_results' for several inputs)")
    parser.add_argument('--version', action='version', version=f'sourmash_taxonomy_parser {__version__}')
    return parser


if __name__ == "__main__":
    parser = build_arg_parser()
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    if not args.sourmash_results or not args.assembly_db or not args.outdir:
        parser.error("the following arguments are required: -s/--sourmash_results, -a/--assembly_db, -o/--outdir")

    setup_logger(log_file=args.log)
    target_taxa = parse_target_taxa(args.target_taxa)

    main(sourmash_files=args.sourmash_results, assembly_dbs=args.assembly_db, target_taxa=target_taxa,
         outdir=args.outdir, exclude_accessions=args.accessions_to_exclude,
         exclude_files=args.exclude_accessions_file, exclude_tolid=args.exclude_tolid,
         exclude_species=args.exclude_species, min_hashes=args.min_hashes, min_containment=args.min_containment,
         margin=args.margin, prefix=args.prefix, contig_sizes=args.contig_sizes, taxdump=args.taxdump,
         div_mapping=args.div_mapping, target_fallback=args.target_fallback,
         min_target_genomes=args.min_target_genomes)
