import pytest
import tempfile
import os
import sys
import pandas as pd          # still used in helper assertions that read CSV outputs
import polars as pl
from collections import defaultdict

# Add bin directory to path to import the script
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import sourmash_taxonomy_parser as stp

# Initialize logger for testing (adds stderr handler; no return value)
stp.setup_logger()

# Sample data for testing
SAMPLE_SOURMASH_CSV = """query_name,match_name,containment,jaccard,intersect_hashes
contig_1,GCA_000001,0.8,0.7,100
contig_1,GCA_000002,0.6,0.5,80
contig_2,GCA_000003,0.9,0.8,120
"""

SAMPLE_ASSEMBLY_CSV = """assembly_accession,taxid,species,genus,family,order,class,phylum,kingdom
GCA_000001,12345,Homo sapiens,Homo,Hominidae,Primates,Mammalia,Chordata,Metazoa
GCA_000002,67890,Mus musculus,Mus,Muridae,Rodentia,Mammalia,Chordata,Metazoa
GCA_000003,11111,Drosophila melanogaster,Drosophila,Drosophilidae,Diptera,Insecta,Arthropoda,Metazoa
"""

@pytest.fixture
def temp_sourmash_file():
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(SAMPLE_SOURMASH_CSV)
        temp_file = f.name
    yield temp_file
    os.unlink(temp_file)

@pytest.fixture
def temp_assembly_file():
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(SAMPLE_ASSEMBLY_CSV)
        temp_file = f.name
    yield temp_file
    os.unlink(temp_file)

@pytest.fixture
def temp_empty_file():
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write("")
        temp_file = f.name
    yield temp_file
    os.unlink(temp_file)

@pytest.fixture
def temp_invalid_sourmash_file():
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write("query_name,match_name,containment,jaccard\ncontig_1,GCA_000001,not_a_number,0.7\n")
        temp_file = f.name
    yield temp_file
    os.unlink(temp_file)

def test_parse_sourmash_results_single_file(temp_sourmash_file):
    """Test parsing a single sourmash results file."""
    result = stp.parse_sourmash_results([temp_sourmash_file])
    assert 'contig_1' in result
    assert 'contig_2' in result
    assert len(result['contig_1']) == 2
    assert len(result['contig_2']) == 1
    assert result['contig_1'][0]['match_name'] == 'GCA_000001'
    assert result['contig_1'][0]['containment'] == 0.8

def test_parse_sourmash_results_multiple_files(temp_sourmash_file):
    """Test parsing multiple sourmash results files and merging."""
    # Create a second file with additional data
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write("query_name,match_name,containment,jaccard,intersect_hashes\ncontig_3,GCA_000004,0.5,0.4,60\n")
        temp_file2 = f.name

    try:
        result = stp.parse_sourmash_results([temp_sourmash_file, temp_file2])
        assert 'contig_1' in result
        assert 'contig_2' in result
        assert 'contig_3' in result
        assert len(result['contig_3']) == 1
    finally:
        os.unlink(temp_file2)

def test_parse_assembly_database(temp_assembly_file):
    """Test parsing assembly database file."""
    result = stp.parse_assembly_database([temp_assembly_file])
    assert isinstance(result, pl.DataFrame)
    assert len(result) == 3
    # assembly_accession is a regular column (not an index) in the polars version
    assert 'GCA_000001' in result.get_column('assembly_accession').to_list()
    homo_row = result.filter(pl.col('assembly_accession') == 'GCA_000001')
    assert homo_row['species'][0] == 'Homo sapiens'

def test_parse_assembly_database_missing_assembly_accession():
    """Assembly DB without 'assembly_accession' column must raise ValueError."""
    csv = "accession,taxid,species\nGCA_1,9606,Homo sapiens\n"
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(csv)
        tmp = f.name
    try:
        with pytest.raises(ValueError, match="assembly_accession"):
            stp.parse_assembly_database([tmp])
    finally:
        os.unlink(tmp)

def test_parse_assembly_database_missing_taxid():
    """Assembly DB without 'taxid' column must raise ValueError."""
    csv = "assembly_accession,species\nGCA_1,Homo sapiens\n"
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(csv)
        tmp = f.name
    try:
        with pytest.raises(ValueError, match="taxid"):
            stp.parse_assembly_database([tmp])
    finally:
        os.unlink(tmp)

def test_get_target_genomes_missing_column_warns(temp_assembly_file, caplog):
    """get_target_genomes warns when requested taxonomic level is absent from DB."""
    import logging
    assembly_df = stp.parse_assembly_database([temp_assembly_file])
    with caplog.at_level(logging.WARNING, logger='sourmash_taxonomy_parser'):
        result = stp.get_target_genomes(assembly_df, {'nonexistent_rank': 'somevalue'})
    assert result == set()
    assert any('nonexistent_rank' in msg for msg in caplog.messages)

def test_generate_summary():
    """Test generating summary from query matches."""
    query_matches = {
        'contig_1': [
            {'match_name': 'GCA_000001', 'containment': 0.8, 'jaccard': 0.7, 'intersect_hashes': 100},
            {'match_name': 'GCA_000002', 'containment': 0.6, 'jaccard': 0.5, 'intersect_hashes': 80}
        ]
    }
    summary = stp.generate_summary(query_matches)
    assert len(summary) == 2
    # Should be sorted by intersect_hashes descending
    assert summary[0][1] == 'GCA_000001'  # higher hashes first
    assert summary[0][2] == 1  # rank
    assert summary[1][2] == 2

def test_write_summary_output(temp_sourmash_file, temp_assembly_file):
    """Test writing summary output."""
    with tempfile.TemporaryDirectory() as temp_dir:
        output_file = os.path.join(temp_dir, 'test_summary.csv')
        query_matches = stp.parse_sourmash_results([temp_sourmash_file])
        summary = stp.generate_summary(query_matches)
        assembly_df = stp.parse_assembly_database([temp_assembly_file])
        target_taxa = {'order': 'primates'}

        stp.write_summary_output(summary, output_file, assembly_df, target_taxa, [temp_sourmash_file], [temp_assembly_file])

        assert os.path.exists(output_file)
        with open(output_file, 'r') as f:
            lines = f.readlines()
            assert len(lines) == 4  # header + 3 data lines
            assert 'header,assembly_accession,taxa,top_n,containment,jaccard,intersect_hashes,is_target' in lines[0]
            # Check no comments (no lines starting with #)

def test_write_non_target_output(temp_sourmash_file, temp_assembly_file):
    """Test writing non-target output."""
    with tempfile.TemporaryDirectory() as temp_dir:
        summary_file = os.path.join(temp_dir, 'summary.csv')
        non_target_file = os.path.join(temp_dir, 'non_target.csv')

        # First create summary file
        query_matches = stp.parse_sourmash_results([temp_sourmash_file])
        summary = stp.generate_summary(query_matches)
        assembly_df = stp.parse_assembly_database([temp_assembly_file])
        target_taxa = {'order': 'primates'}  # Only GCA_000001 matches
        stp.write_summary_output(summary, summary_file, assembly_df, target_taxa, [temp_sourmash_file], [temp_assembly_file])

        # Now create non-target file
        stp.write_non_target_output(summary_file, non_target_file, assembly_df, [temp_sourmash_file], [temp_assembly_file], target_taxa)

        assert os.path.exists(non_target_file)
        with open(non_target_file, 'r') as f:
            lines = f.readlines()
            assert len(lines) >= 2  # header + at least one data line (contig_2 should be non-target)
            assert 'header,assembly_accession,taxa,containment,jaccard,intersect_hashes,species,genus,family,order,class,phylum,kingdom' in lines[0]

def test_main_integration(temp_sourmash_file, temp_assembly_file):
    """Integration test for main function."""
    with tempfile.TemporaryDirectory() as temp_dir:
        stp.main(
            sourmash_files=[temp_sourmash_file],
            assembly_dbs=[temp_assembly_file],
            target_taxa={'order': 'primates'},
            outdir=temp_dir
        )

        # Check if output files are created
        files = os.listdir(temp_dir)
        assert any('summary.csv' in f for f in files)
        assert any('non_target.csv' in f for f in files)

def test_edge_cases_empty_file(temp_empty_file):
    """Test handling of empty files."""
    with pytest.raises(ValueError):
        stp.parse_sourmash_results([temp_empty_file])

def test_edge_cases_invalid_data(temp_invalid_sourmash_file):
    """Test handling of invalid data (non-numeric values)."""
    result = stp.parse_sourmash_results([temp_invalid_sourmash_file])
    # Should skip invalid rows
    assert len(result) == 0  # No valid rows

def test_extract_contig_number():
    """Test extracting contig number from query name."""
    assert stp.extract_contig_number('contig_1') == 1
    assert stp.extract_contig_number('scaffold_42') == 42
    assert stp.extract_contig_number('chr1') == 1
    assert stp.extract_contig_number('no_number') == float('inf')

def test_parse_target_taxa_valid():
    """Test parsing valid target taxa."""
    result = stp.parse_target_taxa(['order:primates', 'family:hominidae'])
    assert result == {'order': 'primates', 'family': 'hominidae'}

def test_parse_target_taxa_invalid():
    """Test parsing invalid target taxa format."""
    result = stp.parse_target_taxa(['orderprimates'])  # No colon
    assert result == {}  # Should be empty due to invalid format

def test_validate_file_paths_exists(temp_sourmash_file):
    """Test validation of existing file paths."""
    # Should not raise exception
    stp.validate_file_paths([temp_sourmash_file], "Test")

def test_validate_file_paths_not_exists():
    """Test validation of non-existing file paths."""
    with pytest.raises(FileNotFoundError):
        stp.validate_file_paths(['/non/existent/file.csv'], "Test")

def test_get_target_genomes_with_matches(temp_assembly_file):
    """Test get_target_genomes with matching taxa."""
    assembly_df = stp.parse_assembly_database([temp_assembly_file])
    target_taxa = {'order': 'primates'}
    result = stp.get_target_genomes(assembly_df, target_taxa)
    assert 'GCA_000001' in result  # Homo sapiens is Primates

def test_get_target_genomes_no_matches(temp_assembly_file):
    """Test get_target_genomes with no matching taxa."""
    assembly_df = stp.parse_assembly_database([temp_assembly_file])
    target_taxa = {'order': 'nonexistent'}
    result = stp.get_target_genomes(assembly_df, target_taxa)
    assert result == set()

def test_parse_sourmash_results_missing_columns():
    """Test parsing sourmash results with missing required columns."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write("query_name,match_name,containment\ncontig_1,GCA_000001,0.8\n")  # Missing jaccard
        temp_file = f.name

    try:
        with pytest.raises(ValueError, match="Missing required columns"):
            stp.parse_sourmash_results([temp_file])
    finally:
        os.unlink(temp_file)

def test_parse_sourmash_results_duplicate_headers():
    """Test parsing sourmash results with duplicate headers (non-numeric rows)."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write("query_name,match_name,containment,jaccard,intersect_hashes\n")
        f.write("query_name,match_name,containment,jaccard,intersect_hashes\n")  # Duplicate header
        f.write("contig_1,GCA_000001,0.8,0.7,100\n")
        temp_file = f.name

    try:
        result = stp.parse_sourmash_results([temp_file])
        assert 'contig_1' in result
        assert len(result['contig_1']) == 1  # Should skip duplicate header row
    finally:
        os.unlink(temp_file)

def test_parse_sourmash_results_with_exclude_accessions(temp_sourmash_file):
    """Test parsing sourmash results with excluded accessions."""
    # Exclude GCA_000001, which should remove one match from contig_1
    result = stp.parse_sourmash_results([temp_sourmash_file], exclude_accessions=['GCA_000001'])
    assert 'contig_1' in result
    assert 'contig_2' in result
    assert len(result['contig_1']) == 1  # Only GCA_000002 remains
    assert len(result['contig_2']) == 1  # GCA_000003 not excluded
    assert result['contig_1'][0]['match_name'] == 'GCA_000002'

def test_parse_sourmash_results_exclude_multiple_accessions(temp_sourmash_file):
    """Test parsing sourmash results with multiple excluded accessions."""
    # Exclude both matches for contig_1
    result = stp.parse_sourmash_results([temp_sourmash_file], exclude_accessions=['GCA_000001', 'GCA_000002'])
    assert 'contig_1' not in result  # No matches left for contig_1
    assert 'contig_2' in result
    assert len(result['contig_2']) == 1

def test_parse_sourmash_results_exclude_nonexistent_accession(temp_sourmash_file):
    """Test parsing sourmash results with excluded accession that doesn't exist."""
    # Exclude a non-existent accession
    result = stp.parse_sourmash_results([temp_sourmash_file], exclude_accessions=['GCA_999999'])
    assert 'contig_1' in result
    assert 'contig_2' in result
    assert len(result['contig_1']) == 2  # No change
    assert len(result['contig_2']) == 1


def test_parse_sourmash_results_exclude_comma_separated_string(temp_sourmash_file):
    """Test that comma-separated accessions in a single string are split correctly.

    This mirrors how ext.args passes the value from Nextflow:
        ext.args = '--accessions_to_exclude GCF_000009045.1,GCF_000013265.1'
    which arrives as a single element list ['GCF_000009045.1,GCF_000013265.1'].
    """
    result = stp.parse_sourmash_results(
        [temp_sourmash_file],
        exclude_accessions=['GCA_000001,GCA_000002']  # single comma-separated string
    )
    assert 'contig_1' not in result   # both matches excluded
    assert 'contig_2' in result       # GCA_000003 untouched


def test_write_summary_output_is_target_values(temp_sourmash_file, temp_assembly_file):
    """Test that is_target column is correctly True/False/None for each row."""
    with tempfile.TemporaryDirectory() as temp_dir:
        output_file = os.path.join(temp_dir, 'summary.csv')
        query_matches = stp.parse_sourmash_results([temp_sourmash_file])
        summary = stp.generate_summary(query_matches)
        assembly_df = stp.parse_assembly_database([temp_assembly_file])
        # Only Primates match → GCA_000001 is target, GCA_000002 and GCA_000003 are not
        target_taxa = {'order': 'primates'}

        stp.write_summary_output(summary, output_file, assembly_df, target_taxa,
                                 [temp_sourmash_file], [temp_assembly_file])

        df = pd.read_csv(output_file)
        # GCA_000001 (Homo sapiens, Primates) → True
        row_target = df[df['assembly_accession'] == 'GCA_000001']
        assert not row_target.empty
        assert row_target.iloc[0]['is_target'] == True

        # GCA_000002 (Rodentia) → False
        row_non = df[df['assembly_accession'] == 'GCA_000002']
        assert not row_non.empty
        assert row_non.iloc[0]['is_target'] == False

        # GCA_000003 (Diptera) → False
        row_non2 = df[df['assembly_accession'] == 'GCA_000003']
        assert not row_non2.empty
        assert row_non2.iloc[0]['is_target'] == False


def test_write_non_target_exact_content(temp_sourmash_file, temp_assembly_file):
    """Test that non_target.csv contains the correct query and lineage info."""
    with tempfile.TemporaryDirectory() as temp_dir:
        summary_file = os.path.join(temp_dir, 'summary.csv')
        non_target_file = os.path.join(temp_dir, 'non_target.csv')

        query_matches = stp.parse_sourmash_results([temp_sourmash_file])
        summary = stp.generate_summary(query_matches)
        assembly_df = stp.parse_assembly_database([temp_assembly_file])
        # order:primates → GCA_000001 is target; contig_1 has a target match so it's excluded from non_target
        # contig_2 top1 = GCA_000003 (Diptera) → non-target
        target_taxa = {'order': 'primates'}
        stp.write_summary_output(summary, summary_file, assembly_df, target_taxa,
                                 [temp_sourmash_file], [temp_assembly_file])
        stp.write_non_target_output(summary_file, non_target_file, assembly_df,
                                    [temp_sourmash_file], [temp_assembly_file], target_taxa)

        df = pd.read_csv(non_target_file)
        # Only contig_2 should appear (contig_1 has a target match)
        assert list(df['header']) == ['contig_2']
        assert df.iloc[0]['assembly_accession'] == 'GCA_000003'
        # Lineage should be populated
        assert df.iloc[0]['order'] == 'Diptera'
        assert df.iloc[0]['phylum'] == 'Arthropoda'


def test_exclude_top1_changes_non_target(temp_assembly_file):
    """When the top1 target match is excluded, the next match becomes top1.

    contig_1 matches: GCA_000001 (Primates, 100 hashes) and GCA_000002 (Rodentia, 80 hashes).
    If GCA_000001 is excluded, GCA_000002 becomes top1 → contig_1 is non-target.
    """
    # Build sourmash results with known hashes so the ranking is deterministic
    sourmash_csv = (
        "query_name,match_name,containment,jaccard,intersect_hashes\n"
        "contig_1,GCA_000001,0.8,0.7,100\n"
        "contig_1,GCA_000002,0.6,0.5,80\n"
    )
    assembly_csv = (
        "assembly_accession,taxid,species,genus,family,order,class,phylum,kingdom\n"
        "GCA_000001,12345,Homo sapiens,Homo,Hominidae,Primates,Mammalia,Chordata,Metazoa\n"
        "GCA_000002,67890,Mus musculus,Mus,Muridae,Rodentia,Mammalia,Chordata,Metazoa\n"
    )

    with tempfile.TemporaryDirectory() as temp_dir:
        sm_file = os.path.join(temp_dir, 'sm.csv')
        db_file = os.path.join(temp_dir, 'db.csv')
        summary_file = os.path.join(temp_dir, 'summary.csv')
        non_target_file = os.path.join(temp_dir, 'non_target.csv')

        with open(sm_file, 'w') as f:
            f.write(sourmash_csv)
        with open(db_file, 'w') as f:
            f.write(assembly_csv)

        # Without exclude: GCA_000001 is top1, is_target=True → contig_1 not in non_target
        query_matches = stp.parse_sourmash_results([sm_file])
        assembly_df = stp.parse_assembly_database([db_file])
        summary = stp.generate_summary(query_matches)
        target_taxa = {'order': 'primates'}
        stp.write_summary_output(summary, summary_file, assembly_df, target_taxa, [sm_file], [db_file])
        stp.write_non_target_output(summary_file, non_target_file, assembly_df, [sm_file], [db_file], target_taxa)

        df_nt = pd.read_csv(non_target_file)
        assert df_nt.empty, "contig_1 should NOT be in non_target when GCA_000001 (target) is top1"

        # With exclude GCA_000001: GCA_000002 becomes top1, is_target=False → contig_1 IS non-target
        summary_excl_file = os.path.join(temp_dir, 'summary_excl.csv')
        non_target_excl_file = os.path.join(temp_dir, 'non_target_excl.csv')

        query_matches_excl = stp.parse_sourmash_results([sm_file], exclude_accessions=['GCA_000001'])
        summary_excl = stp.generate_summary(query_matches_excl)
        stp.write_summary_output(summary_excl, summary_excl_file, assembly_df, target_taxa, [sm_file], [db_file])
        stp.write_non_target_output(summary_excl_file, non_target_excl_file, assembly_df, [sm_file], [db_file], target_taxa)

        df_nt_excl = pd.read_csv(non_target_excl_file)
        assert not df_nt_excl.empty, "contig_1 SHOULD be in non_target after excluding GCA_000001"
        assert df_nt_excl.iloc[0]['header'] == 'contig_1'
        assert df_nt_excl.iloc[0]['assembly_accession'] == 'GCA_000002'


def test_main_integration_with_exclude(temp_sourmash_file, temp_assembly_file):
    """Integration test: main() with exclude_accessions produces fewer summary rows."""
    with tempfile.TemporaryDirectory() as temp_dir:
        stp.main(
            sourmash_files=[temp_sourmash_file],
            assembly_dbs=[temp_assembly_file],
            target_taxa={'order': 'primates'},
            outdir=temp_dir,
            exclude_accessions=['GCA_000001']
        )

        files = os.listdir(temp_dir)
        summary_file = next(f for f in files if 'summary.csv' in f)
        df = pd.read_csv(os.path.join(temp_dir, summary_file))

        # GCA_000001 must not appear anywhere in the output
        assert 'GCA_000001' not in df['assembly_accession'].values
        # contig_1 still present (GCA_000002 remains), contig_2 still present
        assert 'contig_1' in df['header'].values
        assert 'contig_2' in df['header'].values
        # Total rows: contig_1 has 1 match left, contig_2 has 1 → 2 rows
        assert len(df) == 2


# ---------------------------------------------------------------------------
# Tests for unknown-taxid / missing-from-database behaviour
# Rule: if an accession is not in the assembly DB → is_target=True (conservative)
#       so the contig is NEVER written to non_target.csv automatically.
# ---------------------------------------------------------------------------

def _write_pair(tmp, sm_csv, asm_csv, target_taxa):
    """Helper: write CSVs, run full pipeline, return (summary_df, non_target_df or None)."""
    sm_file  = os.path.join(tmp, 'sm.csv')
    db_file  = os.path.join(tmp, 'db.csv')
    with open(sm_file, 'w') as f: f.write(sm_csv)
    with open(db_file, 'w') as f: f.write(asm_csv)

    assembly_df  = stp.parse_assembly_database([db_file])
    query_matches = stp.parse_sourmash_results([sm_file])
    summary      = stp.generate_summary(query_matches)
    summary_file = os.path.join(tmp, 'summary.csv')
    stp.write_summary_output(summary, summary_file, assembly_df, target_taxa, [sm_file], [db_file])

    non_target_file = os.path.join(tmp, 'non_target.csv')
    stp.write_non_target_output(summary_file, non_target_file, assembly_df,
                                [sm_file], [db_file], target_taxa)

    s_df  = pd.read_csv(summary_file)
    nt_df = pd.read_csv(non_target_file)
    return s_df, nt_df


def test_unknown_taxid_treated_as_target_in_summary():
    """Match not present in assembly DB → is_target=True in summary.csv."""
    sm_csv = (
        "query_name,match_name,containment,jaccard,intersect_hashes\n"
        "contig_1,GCA_UNKNOWN,0.9,0.85,900\n"
    )
    asm_csv = (
        "assembly_accession,taxid,species,genus,family,order,class,phylum,kingdom\n"
        "GCA_000001,9606,Homo sapiens,Homo,Hominidae,Primates,Mammalia,Chordata,Animalia\n"
    )
    with tempfile.TemporaryDirectory() as tmp:
        s_df, _ = _write_pair(tmp, sm_csv, asm_csv, {'order': 'primates'})
        row = s_df[s_df['assembly_accession'] == 'GCA_UNKNOWN']
        assert not row.empty
        # Unknown accession must be flagged as target (conservative — never auto-remove)
        assert row.iloc[0]['is_target'] == True


def test_unknown_taxid_contig_not_in_non_target():
    """Contig whose only match is unknown-taxid must NOT appear in non_target.csv."""
    sm_csv = (
        "query_name,match_name,containment,jaccard,intersect_hashes\n"
        "contig_1,GCA_UNKNOWN,0.9,0.85,900\n"
    )
    asm_csv = (
        "assembly_accession,taxid,species,genus,family,order,class,phylum,kingdom\n"
        "GCA_000001,9606,Homo sapiens,Homo,Hominidae,Primates,Mammalia,Chordata,Animalia\n"
    )
    with tempfile.TemporaryDirectory() as tmp:
        _, nt_df = _write_pair(tmp, sm_csv, asm_csv, {'order': 'primates'})
        assert 'contig_1' not in nt_df['header'].values, (
            "contig_1 matched only an unknown accession → should be treated as target, "
            "NOT written to non_target.csv"
        )


def test_mixed_known_unknown_contig_not_in_non_target():
    """Contig with one known non-target match AND one unknown match → treated as target.

    Because the unknown match gets is_target=True, the contig has at least one
    target match and must NOT appear in non_target.csv.
    """
    sm_csv = (
        "query_name,match_name,containment,jaccard,intersect_hashes\n"
        "contig_1,GCA_000002,0.6,0.5,60\n"   # known, Rodentia → non-target vs primates
        "contig_1,GCA_UNKNOWN,0.8,0.7,80\n"  # unknown → treated as target
    )
    asm_csv = (
        "assembly_accession,taxid,species,genus,family,order,class,phylum,kingdom\n"
        "GCA_000001,9606,Homo sapiens,Homo,Hominidae,Primates,Mammalia,Chordata,Animalia\n"
        "GCA_000002,10090,Mus musculus,Mus,Muridae,Rodentia,Mammalia,Chordata,Animalia\n"
    )
    with tempfile.TemporaryDirectory() as tmp:
        _, nt_df = _write_pair(tmp, sm_csv, asm_csv, {'order': 'primates'})
        assert 'contig_1' not in nt_df['header'].values, (
            "contig_1 has an unknown match (treated as target) → must NOT be in non_target"
        )


def test_all_nontarget_known_contig_in_non_target():
    """Contig with exclusively known non-target matches → appears in non_target.csv."""
    sm_csv = (
        "query_name,match_name,containment,jaccard,intersect_hashes\n"
        "contig_1,GCA_000002,0.7,0.6,70\n"
        "contig_1,GCA_000003,0.5,0.4,50\n"
    )
    asm_csv = (
        "assembly_accession,taxid,species,genus,family,order,class,phylum,kingdom\n"
        "GCA_000002,10090,Mus musculus,Mus,Muridae,Rodentia,Mammalia,Chordata,Animalia\n"
        "GCA_000003,7227,Drosophila melanogaster,Drosophila,Drosophilidae,Diptera,Insecta,Arthropoda,Animalia\n"
    )
    with tempfile.TemporaryDirectory() as tmp:
        _, nt_df = _write_pair(tmp, sm_csv, asm_csv, {'order': 'primates'})
        assert 'contig_1' in nt_df['header'].values, (
            "contig_1 has only known non-target matches → must appear in non_target.csv"
        )


def test_single_target_match_prevents_non_target_entry():
    """Even one target match among many non-target ones → contig excluded from non_target.csv.

    This is the strict rule: ALL matches must be non-target.
    """
    sm_csv = (
        "query_name,match_name,containment,jaccard,intersect_hashes\n"
        "contig_1,GCA_000001,0.9,0.85,900\n"   # target (Primates)
        "contig_1,GCA_000002,0.7,0.6,70\n"     # non-target (Rodentia)
        "contig_1,GCA_000003,0.5,0.4,50\n"     # non-target (Diptera)
    )
    asm_csv = (
        "assembly_accession,taxid,species,genus,family,order,class,phylum,kingdom\n"
        "GCA_000001,9606,Homo sapiens,Homo,Hominidae,Primates,Mammalia,Chordata,Animalia\n"
        "GCA_000002,10090,Mus musculus,Mus,Muridae,Rodentia,Mammalia,Chordata,Animalia\n"
        "GCA_000003,7227,Drosophila melanogaster,Drosophila,Drosophilidae,Diptera,Insecta,Arthropoda,Animalia\n"
    )
    with tempfile.TemporaryDirectory() as tmp:
        _, nt_df = _write_pair(tmp, sm_csv, asm_csv, {'order': 'primates'})
        assert 'contig_1' not in nt_df['header'].values, (
            "contig_1 has at least one target match → must NOT appear in non_target.csv "
            "(strict rule: ALL matches must be non-target)"
        )


def test_non_target_output_top1_row_written():
    """The row written to non_target.csv is the top1 match (highest intersect_hashes)."""
    sm_csv = (
        "query_name,match_name,containment,jaccard,intersect_hashes\n"
        "contig_1,GCA_000002,0.5,0.4,50\n"   # rank 2 (fewer hashes)
        "contig_1,GCA_000003,0.7,0.6,200\n"  # rank 1 (more hashes) → this should appear
    )
    asm_csv = (
        "assembly_accession,taxid,species,genus,family,order,class,phylum,kingdom\n"
        "GCA_000002,10090,Mus musculus,Mus,Muridae,Rodentia,Mammalia,Chordata,Animalia\n"
        "GCA_000003,7227,Drosophila melanogaster,Drosophila,Drosophilidae,Diptera,Insecta,Arthropoda,Animalia\n"
    )
    with tempfile.TemporaryDirectory() as tmp:
        _, nt_df = _write_pair(tmp, sm_csv, asm_csv, {'order': 'primates'})
        assert 'contig_1' in nt_df['header'].values
        row = nt_df[nt_df['header'] == 'contig_1'].iloc[0]
        assert row['assembly_accession'] == 'GCA_000003', (
            "non_target.csv should contain the top1 match (most intersect_hashes)"
        )
