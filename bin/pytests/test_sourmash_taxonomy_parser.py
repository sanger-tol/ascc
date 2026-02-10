import pytest
import tempfile
import os
import sys
import pandas as pd
from collections import defaultdict

# Add bin directory to path to import the script
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
import sourmash_taxonomy_parser as stp

# Initialize logger for testing
stp.logger = stp.setup_logger()

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
    assert isinstance(result, pd.DataFrame)
    assert len(result) == 3
    assert 'GCA_000001' in result.index
    assert result.loc['GCA_000001', 'species'] == 'Homo sapiens'

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
    result = stp.parse_target_taxa(['order:primates', 'family:hominidae'], stp.logger)
    assert result == {'order': 'primates', 'family': 'hominidae'}

def test_parse_target_taxa_invalid():
    """Test parsing invalid target taxa format."""
    result = stp.parse_target_taxa(['orderprimates'], stp.logger)  # No colon
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