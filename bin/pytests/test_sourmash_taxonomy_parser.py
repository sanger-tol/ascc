import math
import os
import sys
import tempfile

import pandas as pd          # used in helper assertions that read CSV outputs
import polars as pl
import pytest

# Add bin directory to path to import the script
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import sourmash_taxonomy_parser as stp  # noqa: E402

# Initialize logger for testing (adds stderr handler; no return value)
stp.setup_logger()

# Sample data for testing (no ksize/scaled columns: old/hand-made format, scaled defaults to 200)
SAMPLE_SOURMASH_CSV = """query_name,match_name,containment,jaccard,intersect_hashes
contig_1,GCA_000001,0.8,0.7,100
contig_1,GCA_000002,0.6,0.5,80
contig_2,GCA_000003,0.9,0.8,120
"""

SAMPLE_ASSEMBLY_CSV = """assembly_accession,assembly_name,taxid,species,genus,family,order,class,phylum,kingdom
GCA_000001,hsHomSapi1.hap1.1,12345,Homo sapiens,Homo,Hominidae,Primates,Mammalia,Chordata,Metazoa
GCA_000002,ASM_mouse,67890,Mus musculus,Mus,Muridae,Rodentia,Mammalia,Chordata,Metazoa
GCA_000003,idDroMela1.1,11111,Drosophila melanogaster,Drosophila,Drosophilidae,Diptera,Insecta,Arthropoda,Metazoa
"""

ASM_HEADER = "assembly_accession,taxid,species,genus,family,order,class,phylum,kingdom\n"
HOMO = "GCA_000001,9606,Homo sapiens,Homo,Hominidae,Primates,Mammalia,Chordata,Animalia\n"
MUS = "GCA_000002,10090,Mus musculus,Mus,Muridae,Rodentia,Mammalia,Chordata,Animalia\n"
DRO = "GCA_000003,7227,Drosophila melanogaster,Drosophila,Drosophilidae,Diptera,Insecta,Arthropoda,Animalia\n"


def _tmp(content, suffix='.csv'):
    with tempfile.NamedTemporaryFile(mode='w', suffix=suffix, delete=False) as f:
        f.write(content)
        return f.name


@pytest.fixture
def temp_sourmash_file():
    p = _tmp(SAMPLE_SOURMASH_CSV)
    yield p
    os.unlink(p)


@pytest.fixture
def temp_assembly_file():
    p = _tmp(SAMPLE_ASSEMBLY_CSV)
    yield p
    os.unlink(p)


def run(tmp, sm_csv, asm_csv, target_taxa, **kw):
    """Write CSVs, run main(), return (summary_df, non_target_df, classification_df)."""
    sm_file = os.path.join(tmp, 'sm.csv')
    db_file = os.path.join(tmp, 'db.csv')
    with open(sm_file, 'w') as f:
        f.write(sm_csv)
    with open(db_file, 'w') as f:
        f.write(asm_csv)
    outdir = os.path.join(tmp, 'out')
    stp.main(sourmash_files=[sm_file], assembly_dbs=[db_file], target_taxa=target_taxa, outdir=outdir, **kw)
    files = os.listdir(outdir)
    # keep_default_na=False: 'None' / '' stay strings as downstream consumers see them
    read = lambda sfx: pd.read_csv(os.path.join(outdir, next(f for f in files if f.endswith(sfx))),  # noqa: E731
                                   keep_default_na=False)
    nt = read('.non_target.csv') if any(f.endswith('.non_target.csv') for f in files) else None
    return read('.summary.csv'), nt, read('.classification.csv')


def cls(c_df, contig):
    """Readable label from the FCS-GX-style classification row.

    EXCLUDE -> 'non_target', REVIEW -> 'review', NA split by evidence: 'no_hit' (no hits left),
    'target' (reported hit is a target), 'unclassified' (only weak non-target hits)."""
    row = c_df.set_index('scaff').loc[contig]
    if row['sourmash_action'] == 'EXCLUDE':
        return 'non_target'
    if row['sourmash_action'] == 'REVIEW':
        return 'review'
    assert row['sourmash_action'] == 'NA'
    if int(row['sourmash_n_hits']) == 0:
        return 'no_hit'
    return 'target' if str(row['sourmash_top_is_target']).lower() == 'true' else 'unclassified'


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def test_parse_sourmash_results_single_file(temp_sourmash_file):
    hits, queries = stp.parse_sourmash_results([temp_sourmash_file])
    assert isinstance(hits, pl.DataFrame)
    assert sorted(queries) == ['contig_1', 'contig_2']
    assert hits.filter(pl.col('query_name') == 'contig_1').height == 2
    # scaled defaults to 200 when the column is absent
    assert hits.get_column('scaled').unique().to_list() == [200]


def test_parse_sourmash_results_multiple_files(temp_sourmash_file):
    f2 = _tmp("query_name,match_name,containment,jaccard,intersect_hashes,ksize,scaled\ncontig_3,GCA_000004 some virus,0.5,0.4,60,21,50\n")
    try:
        hits, queries = stp.parse_sourmash_results([temp_sourmash_file, f2])
        assert 'contig_3' in queries
        row = hits.filter(pl.col('query_name') == 'contig_3').row(0, named=True)
        assert row['match_name'] == 'GCA_000004'   # description stripped
        assert row['ksize'] == 21 and row['scaled'] == 50
    finally:
        os.unlink(f2)


def test_db_label_from_multisearch_file_name():
    assert stp.db_label('x/icCatNigr1_HAP1_vs_gtdb_v226_k31_s200.multisearch_results.csv') == 'gtdb_v226_k31_s200'
    assert stp.db_label('merged.csv') == 'merged.csv'


def test_parse_assembly_database(temp_assembly_file):
    result = stp.parse_assembly_database([temp_assembly_file])
    assert isinstance(result, pl.DataFrame)
    assert len(result) == 3
    assert result.filter(pl.col('assembly_accession') == 'GCA_000001')['species'][0] == 'Homo sapiens'


def test_parse_assembly_database_missing_assembly_accession():
    p = _tmp("accession,taxid,species\nGCA_1,9606,Homo sapiens\n")
    try:
        with pytest.raises(ValueError, match="assembly_accession"):
            stp.parse_assembly_database([p])
    finally:
        os.unlink(p)


def test_parse_assembly_database_missing_taxid():
    p = _tmp("assembly_accession,species\nGCA_1,Homo sapiens\n")
    try:
        with pytest.raises(ValueError, match="taxid"):
            stp.parse_assembly_database([p])
    finally:
        os.unlink(p)


def test_edge_cases_empty_file():
    p = _tmp("")
    try:
        with pytest.raises(ValueError):
            stp.parse_sourmash_results([p])
    finally:
        os.unlink(p)


def test_edge_cases_invalid_data():
    p = _tmp("query_name,match_name,containment,jaccard\ncontig_1,GCA_000001,not_a_number,0.7\n")
    try:
        hits, _ = stp.parse_sourmash_results([p])
        assert hits.height == 0
    finally:
        os.unlink(p)


def test_parse_sourmash_results_missing_columns():
    p = _tmp("query_name,match_name,containment\ncontig_1,GCA_000001,0.8\n")
    try:
        with pytest.raises(ValueError, match="Missing required columns"):
            stp.parse_sourmash_results([p])
    finally:
        os.unlink(p)


def test_parse_sourmash_results_duplicate_headers():
    p = _tmp("query_name,match_name,containment,jaccard,intersect_hashes\n"
             "query_name,match_name,containment,jaccard,intersect_hashes\n"
             "contig_1,GCA_000001,0.8,0.7,100\n")
    try:
        hits, _ = stp.parse_sourmash_results([p])
        assert hits.height == 1
    finally:
        os.unlink(p)


def test_extract_contig_number():
    assert stp.extract_contig_number('contig_1') == 1
    assert stp.extract_contig_number('scaffold_42') == 42
    assert stp.extract_contig_number('no_number') == float('inf')


def test_parse_target_taxa():
    assert stp.parse_target_taxa(['order:primates', 'family:hominidae']) == {'order': 'primates', 'family': 'hominidae'}
    assert stp.parse_target_taxa(['orderprimates']) == {}


def test_validate_file_paths():
    with pytest.raises(FileNotFoundError):
        stp.validate_file_paths(['/non/existent/file.csv'], "Test")


def test_get_target_genomes(temp_assembly_file, caplog):
    import logging
    df = stp.parse_assembly_database([temp_assembly_file])
    assert 'GCA_000001' in stp.get_target_genomes(df, {'order': 'primates'})
    assert stp.get_target_genomes(df, {'order': 'nonexistent'}) is None   # sentinel: taxon absent
    assert stp.get_target_genomes(df, {}) == set()
    with caplog.at_level(logging.WARNING, logger='sourmash_taxonomy_parser'):
        stp.get_target_genomes(df, {'nonexistent_rank': 'x'})
    assert any('nonexistent_rank' in m for m in caplog.messages)


# ---------------------------------------------------------------------------
# Exclusions
# ---------------------------------------------------------------------------

def test_exclude_accessions_cli_list_and_comma(temp_sourmash_file):
    ex = stp.build_exclude_set(['GCA_000001,GCA_000002'])
    hits, queries = stp.parse_sourmash_results([temp_sourmash_file], exclude_set=ex)
    assert 'contig_1' not in hits.get_column('query_name').to_list()
    assert 'contig_1' in queries          # still known, reported as no_hit later


def test_exclude_accessions_backward_compatible_argument(temp_sourmash_file):
    hits, _ = stp.parse_sourmash_results([temp_sourmash_file], exclude_accessions=['GCA_000001'])
    assert 'GCA_000001' not in hits.get_column('match_name').to_list()
    assert hits.height == 2


def test_exclude_unversioned_matches_any_version():
    p = _tmp("query_name,match_name,containment,jaccard,intersect_hashes\n"
             "c1,GCA_000001.1,0.8,0.7,100\nc1,GCA_000001.2,0.8,0.7,100\nc1,GCA_000011.1,0.8,0.7,100\n")
    try:
        hits, _ = stp.parse_sourmash_results([p], exclude_set={'GCA_000001'})
        assert hits.get_column('match_name').to_list() == ['GCA_000011.1']
        hits, _ = stp.parse_sourmash_results([p], exclude_set={'GCA_000001.2'})
        assert sorted(hits.get_column('match_name').to_list()) == ['GCA_000001.1', 'GCA_000011.1']
    finally:
        os.unlink(p)


def test_read_exclude_file_tsv_with_column():
    p = _tmp("assembly_accession\tdb\ttax_name\treason\nGCA_1.1\tvertebrate\tX y\tncbi:contaminated\nGCA_2.1\tviruses\tZ\tmanual\n", '.tsv')
    try:
        assert stp.read_exclude_file(p) == {'GCA_1.1', 'GCA_2.1'}
    finally:
        os.unlink(p)


def test_read_exclude_file_csv_with_column_not_first():
    p = _tmp("reason,assembly_accession\nmanual,GCA_3.1\n")
    try:
        assert stp.read_exclude_file(p) == {'GCA_3.1'}
    finally:
        os.unlink(p)


def test_read_exclude_file_plain_list():
    p = _tmp("# comment\nGCA_1.1\n\nGCA_2.1  some note\n", '.txt')
    try:
        assert stp.read_exclude_file(p) == {'GCA_1.1', 'GCA_2.1'}
    finally:
        os.unlink(p)


def test_tolid_accessions(temp_assembly_file):
    df = stp.parse_assembly_database([temp_assembly_file])
    assert stp.tolid_accessions(df, 'hsHomSapi1') == {'GCA_000001'}
    assert stp.tolid_accessions(df, 'hsHomSapi') == set()       # not a ToLID → skipped
    assert stp.tolid_accessions(df, 'hsHomSapi10') == set()     # other specimen, no prefix match
    assert stp.tolid_accessions(df, 'asccTinyTest_V2') == set()


def test_exclude_tolid_in_main():
    sm = ("query_name,match_name,containment,jaccard,intersect_hashes\n"
          "contig_1,GCA_000001,0.9,0.8,500\n"
          "contig_1,GCA_000003,0.5,0.4,300\n")
    asm = ("assembly_accession,assembly_name,taxid,species,genus,family,order,class,phylum,kingdom\n"
           "GCA_000001,hsHomSapi1.hap1.1,9606,Homo sapiens,Homo,Hominidae,Primates,Mammalia,Chordata,Animalia\n"
           "GCA_000003,x,7227,Drosophila melanogaster,Drosophila,Drosophilidae,Diptera,Insecta,Arthropoda,Animalia\n")
    with tempfile.TemporaryDirectory() as tmp:
        s_df, nt_df, c_df = run(tmp, sm, asm, {'order': 'primates'}, exclude_tolid='hsHomSapi1')
        assert 'GCA_000001' not in s_df['assembly_accession'].values
        assert cls(c_df, 'contig_1') == 'non_target'


def test_excluded_only_contig_is_no_hit():
    sm = ("query_name,match_name,containment,jaccard,intersect_hashes\n"
          "contig_1,GCA_000001,0.9,0.8,500\n"
          "contig_2,GCA_000003,0.5,0.4,300\n")
    with tempfile.TemporaryDirectory() as tmp:
        _, nt_df, c_df = run(tmp, sm, ASM_HEADER + HOMO + DRO, {'order': 'primates'}, exclude_accessions=['GCA_000001'])
        assert cls(c_df, 'contig_1') == 'no_hit'
        assert c_df.set_index('scaff').loc['contig_1', 'sourmash_n_hits'] == 0
        assert 'contig_1' not in nt_df['header'].values


# ---------------------------------------------------------------------------
# Scoring / ranking
# ---------------------------------------------------------------------------

def test_ranking_by_containment_not_hashes():
    """A big virus hit (many s50 hashes, low containment) must not outrank a better euk hit."""
    sm = ("query_name,match_name,containment,jaccard,intersect_hashes,ksize,scaled\n"
          "contig_1,GCA_000003,0.0015,0.0001,600,21,50\n"
          "contig_1,GCA_000001,0.02,0.01,16,51,200\n")
    with tempfile.TemporaryDirectory() as tmp:
        s_df, _, c_df = run(tmp, sm, ASM_HEADER + HOMO + DRO, {'order': 'primates'})
        top = s_df[s_df['top_n'] == 1].iloc[0]
        assert top['assembly_accession'] == 'GCA_000001'
        # norm_hashes = hashes * scaled / 200
        virus = s_df[s_df['assembly_accession'] == 'GCA_000003'].iloc[0]
        assert virus['norm_hashes'] == pytest.approx(150.0)
        assert virus['ksize'] == 21 and virus['scaled'] == 50
        assert cls(c_df, 'contig_1') == 'target'


def test_summary_columns_backward_compatible(temp_sourmash_file, temp_assembly_file):
    with tempfile.TemporaryDirectory() as tmp:
        stp.main(sourmash_files=[temp_sourmash_file], assembly_dbs=[temp_assembly_file],
                 target_taxa={'order': 'primates'}, outdir=tmp)
        files = os.listdir(tmp)
        summary = next(f for f in files if f.endswith('summary.csv'))
        with open(os.path.join(tmp, summary)) as fh:
            header = fh.readline().strip().split(',')
        assert header[:8] == ["header", "assembly_accession", "taxa", "top_n", "containment", "jaccard", "intersect_hashes", "is_target"]
        assert {'db', 'ksize', 'scaled', 'norm_hashes', 'passes_thresholds'} <= set(header)
        nt = next(f for f in files if f.endswith('non_target.csv'))
        with open(os.path.join(tmp, nt)) as fh:
            header = fh.readline().strip().split(',')
        assert header[:13] == ["header", "assembly_accession", "taxa", "containment", "jaccard", "intersect_hashes",
                               "species", "genus", "family", "order", "class", "phylum", "kingdom"]
        assert any(f.endswith('classification.csv') for f in files)


def test_prefix_output_name(temp_sourmash_file, temp_assembly_file):
    with tempfile.TemporaryDirectory() as tmp:
        stp.main(sourmash_files=[temp_sourmash_file], assembly_dbs=[temp_assembly_file],
                 target_taxa={'order': 'primates'}, outdir=tmp, prefix='sourmash_results')
        assert sorted(os.listdir(tmp)) == ['sourmash_results_order_primates.classification.csv',
                                          'sourmash_results_order_primates.non_target.csv',
                                          'sourmash_results_order_primates.summary.csv']


def test_is_target_values(temp_sourmash_file, temp_assembly_file):
    with tempfile.TemporaryDirectory() as tmp:
        s_df, nt_df, c_df = run(tmp, SAMPLE_SOURMASH_CSV, SAMPLE_ASSEMBLY_CSV, {'order': 'primates'})
        by_acc = s_df.set_index('assembly_accession')['is_target']
        assert by_acc['GCA_000001'] == True    # noqa: E712
        assert by_acc['GCA_000002'] == False   # noqa: E712
        assert by_acc['GCA_000003'] == False   # noqa: E712
        # contig_1: target 0.8 beats Rodentia 0.6; contig_2: Diptera, 120 hashes → non-target
        assert cls(c_df, 'contig_1') == 'target'
        assert list(nt_df['header']) == ['contig_2']
        assert nt_df.iloc[0]['order'] == 'Diptera'
        assert nt_df.iloc[0]['phylum'] == 'Arthropoda'


# ---------------------------------------------------------------------------
# Thresholds and margin
# ---------------------------------------------------------------------------

def test_single_hash_non_target_is_unclassified():
    """One shared hash is not evidence: the contig is unclassified (kept), not non-target."""
    sm = ("query_name,match_name,containment,jaccard,intersect_hashes,ksize,scaled\n"
          "contig_1,GCA_000003,0.01,0.001,1,31,200\n")
    with tempfile.TemporaryDirectory() as tmp:
        s_df, nt_df, c_df = run(tmp, sm, ASM_HEADER + HOMO + DRO, {'order': 'primates'}, min_hashes=5)
        assert cls(c_df, 'contig_1') == 'unclassified'
        assert nt_df.empty
        assert s_df.iloc[0]['passes_thresholds'] == False   # noqa: E712


def test_min_hashes_normalised_to_scaled_200():
    """20 hashes at scaled=50 = 5 hashes at scaled=200."""
    sm = ("query_name,match_name,containment,jaccard,intersect_hashes,ksize,scaled\n"
          "contig_1,GCA_000003,0.01,0.001,20,21,50\n")
    with tempfile.TemporaryDirectory() as tmp:
        _, _, c_df = run(tmp, sm, ASM_HEADER + HOMO + DRO, {'order': 'primates'}, min_hashes=5)
        assert cls(c_df, 'contig_1') == 'non_target'
    with tempfile.TemporaryDirectory() as tmp:
        _, _, c_df = run(tmp, sm, ASM_HEADER + HOMO + DRO, {'order': 'primates'}, min_hashes=6)
        assert cls(c_df, 'contig_1') == 'unclassified'


def test_min_containment():
    sm = ("query_name,match_name,containment,jaccard,intersect_hashes\n"
          "contig_1,GCA_000003,0.004,0.001,50\n")
    with tempfile.TemporaryDirectory() as tmp:
        _, _, c_df = run(tmp, sm, ASM_HEADER + HOMO + DRO, {'order': 'primates'}, min_containment=0.005)
        assert cls(c_df, 'contig_1') == 'unclassified'


def test_default_min_containment_protects_large_scaffold():
    """A chromosome-scale scaffold sharing a few dozen hashes (containment < 0.1%) is not called non-target."""
    sm = ("query_name,match_name,containment,jaccard,intersect_hashes,ksize,scaled\n"
          "scaffold_1,GCA_000003,0.00025,0.00001,88,21,50\n"      # 22 s200-equivalent hashes
          "scaffold_1,GCA_000001,0.00012,0.00001,10,51,200\n")
    with tempfile.TemporaryDirectory() as tmp:
        _, _, c_df = run(tmp, sm, ASM_HEADER + HOMO + DRO, {'order': 'primates'})
        assert cls(c_df, 'scaffold_1') == 'target'
    with tempfile.TemporaryDirectory() as tmp:
        _, _, c_df = run(tmp, sm, ASM_HEADER + HOMO + DRO, {'order': 'primates'}, min_containment=0)
        assert cls(c_df, 'scaffold_1') == 'non_target'


def test_margin():
    """EXCLUDE needs score > margin * best target score; beating the target by less is REVIEW."""
    sm = ("query_name,match_name,containment,jaccard,intersect_hashes\n"
          "contig_1,GCA_000001,0.10,0.05,20\n"     # target
          "contig_1,GCA_000003,0.15,0.08,30\n")    # non-target, 1.5x the target
    asm = ASM_HEADER + HOMO + DRO
    with tempfile.TemporaryDirectory() as tmp:
        _, nt_df, c_df = run(tmp, sm, asm, {'order': 'primates'}, margin=1)
        assert cls(c_df, 'contig_1') == 'non_target'
        row = nt_df.iloc[0]
        assert row['assembly_accession'] == 'GCA_000003'
        assert row['best_target_accession'] == 'GCA_000001'
    with tempfile.TemporaryDirectory() as tmp:
        _, nt_df, c_df = run(tmp, sm, asm, {'order': 'primates'}, margin=2)
        assert cls(c_df, 'contig_1') == 'review'
        assert nt_df.empty   # REVIEW is not written to non_target.csv
        row = c_df.set_index('scaff').loc['contig_1']
        assert row['sourmash_top_accession'] == 'GCA_000003'           # the non-target evidence is reported
        assert row['sourmash_best_target_accession'] == 'GCA_000001'
    with tempfile.TemporaryDirectory() as tmp:
        _, _, c_df = run(tmp, sm, asm, {'order': 'primates'}, margin=5)
        assert cls(c_df, 'contig_1') == 'review'
    sm_weaker = sm.replace("0.15,0.08,30", "0.08,0.04,30")               # non-target below the target
    with tempfile.TemporaryDirectory() as tmp:
        _, _, c_df = run(tmp, sm_weaker, asm, {'order': 'primates'})
        assert cls(c_df, 'contig_1') == 'target'


def test_weak_non_target_does_not_beat_target():
    """Best non-target hit above the target by score but without evidence → target."""
    sm = ("query_name,match_name,containment,jaccard,intersect_hashes\n"
          "contig_1,GCA_000003,0.05,0.001,1\n"
          "contig_1,GCA_000001,0.01,0.001,10\n")
    with tempfile.TemporaryDirectory() as tmp:
        s_df, _, c_df = run(tmp, sm, ASM_HEADER + HOMO + DRO, {'order': 'primates'})
        assert s_df[s_df['top_n'] == 1].iloc[0]['assembly_accession'] == 'GCA_000003'
        assert cls(c_df, 'contig_1') == 'target'


def test_non_target_row_is_best_evidence_hit():
    """non_target.csv reports the best non-target hit that passes thresholds, not a weaker top-1."""
    sm = ("query_name,match_name,containment,jaccard,intersect_hashes\n"
          "contig_1,GCA_000002,0.5,0.4,2\n"      # top-1 by containment but only 2 hashes
          "contig_1,GCA_000003,0.3,0.2,200\n")
    with tempfile.TemporaryDirectory() as tmp:
        _, nt_df, _ = run(tmp, sm, ASM_HEADER + HOMO + MUS + DRO, {'order': 'primates'}, min_hashes=5)
        assert nt_df.iloc[0]['assembly_accession'] == 'GCA_000003'


def test_legacy_rule_reproduced():
    """--min_hashes 0 --min_containment 0 --margin inf = parser 1.x EXCLUDE set: any target hit keeps the
    contig (a stronger non-target hit then only makes it REVIEW)."""
    sm = ("query_name,match_name,containment,jaccard,intersect_hashes\n"
          "contig_1,GCA_000001,0.001,0.001,1\n"
          "contig_1,GCA_000003,0.9,0.8,900\n"
          "contig_2,GCA_000003,0.001,0.001,1\n")
    with tempfile.TemporaryDirectory() as tmp:
        _, _, c_df = run(tmp, sm, ASM_HEADER + HOMO + DRO, {'order': 'primates'},
                         min_hashes=0, min_containment=0, margin=math.inf)
        assert cls(c_df, 'contig_1') == 'review'
        assert cls(c_df, 'contig_2') == 'non_target'
    with tempfile.TemporaryDirectory() as tmp:
        _, _, c_df = run(tmp, sm, ASM_HEADER + HOMO + DRO, {'order': 'primates'})
        assert cls(c_df, 'contig_1') == 'non_target'
        assert cls(c_df, 'contig_2') == 'unclassified'


# ---------------------------------------------------------------------------
# Conservative rules (unknown taxid, taxon absent from DB, no target taxa)
# ---------------------------------------------------------------------------

def test_unknown_taxid_treated_as_target():
    sm = "query_name,match_name,containment,jaccard,intersect_hashes\ncontig_1,GCA_UNKNOWN,0.9,0.85,900\n"
    with tempfile.TemporaryDirectory() as tmp:
        s_df, nt_df, c_df = run(tmp, sm, ASM_HEADER + HOMO, {'order': 'primates'})
        assert s_df.iloc[0]['is_target'] == True   # noqa: E712
        assert cls(c_df, 'contig_1') == 'target'
        assert nt_df.empty


def test_unknown_match_protects_contig_with_weaker_non_target():
    sm = ("query_name,match_name,containment,jaccard,intersect_hashes\n"
          "contig_1,GCA_000002,0.6,0.5,60\n"
          "contig_1,GCA_UNKNOWN,0.8,0.7,80\n")
    with tempfile.TemporaryDirectory() as tmp:
        _, nt_df, _ = run(tmp, sm, ASM_HEADER + HOMO + MUS, {'order': 'primates'})
        assert 'contig_1' not in nt_df['header'].values


def test_all_non_target_known_contig_in_non_target():
    sm = ("query_name,match_name,containment,jaccard,intersect_hashes\n"
          "contig_1,GCA_000002,0.7,0.6,70\n"
          "contig_1,GCA_000003,0.5,0.4,50\n")
    with tempfile.TemporaryDirectory() as tmp:
        _, nt_df, _ = run(tmp, sm, ASM_HEADER + HOMO + MUS + DRO, {'order': 'primates'})
        assert list(nt_df['header']) == ['contig_1']
        assert nt_df.iloc[0]['assembly_accession'] == 'GCA_000002'


def test_target_taxon_absent_from_db_keeps_everything():
    with tempfile.TemporaryDirectory() as tmp:
        s_df, nt_df, c_df = run(tmp, SAMPLE_SOURMASH_CSV, SAMPLE_ASSEMBLY_CSV, {'phylum': 'bryophyta'})
        assert (s_df['is_target'] == True).all()   # noqa: E712
        assert nt_df.empty
        assert {cls(c_df, c) for c in c_df['scaff']} == {'target'}


def test_no_target_taxa_writes_summary_and_classification_only(temp_sourmash_file, temp_assembly_file):
    with tempfile.TemporaryDirectory() as tmp:
        s_df, nt_df, c_df = run(tmp, SAMPLE_SOURMASH_CSV, SAMPLE_ASSEMBLY_CSV, {})
        assert nt_df is None
        assert (s_df['is_target'].astype(str) == 'None').all()
        assert set(c_df['sourmash_action']) == {'NA'}


def test_species_falls_back_to_tax_name():
    sm = "query_name,match_name,containment,jaccard,intersect_hashes\ncontig_1,GCA_V,0.5,0.4,50\n"
    asm = ("assembly_accession,taxid,tax_name,species,genus,family,order,class,phylum,kingdom\n"
           "GCA_V,340016,uncultured virus,,,,,,,\n")
    with tempfile.TemporaryDirectory() as tmp:
        _, nt_df, _ = run(tmp, sm, asm + HOMO.replace('GCA_000001,9606,', 'GCA_000001,9606,Homo sapiens,'),
                          {'order': 'primates'})
        assert nt_df.iloc[0]['species'] == 'uncultured virus'
        assert nt_df.iloc[0]['genus'] == 'None'


# ---------------------------------------------------------------------------
# FCS-GX-style classification.csv: columns, all contigs, gx-div, own-species exclusion
# ---------------------------------------------------------------------------

def test_classification_columns_mirror_parsed_fcsgx():
    with tempfile.TemporaryDirectory() as tmp:
        _, _, c_df = run(tmp, SAMPLE_SOURMASH_CSV, SAMPLE_ASSEMBLY_CSV, {'order': 'primates'})
        assert list(c_df.columns) == stp.CLASSIFICATION_COLUMNS
        assert list(c_df.columns[:5]) == ['scaff', 'sourmash_top_tax_name', 'sourmash_top_taxid', 'sourmash_div', 'sourmash_action']
        assert set(c_df['sourmash_action']) <= {'EXCLUDE', 'REVIEW', 'NA'}
        row = c_df.set_index('scaff').loc['contig_2']
        assert row['sourmash_action'] == 'EXCLUDE'
        assert row['sourmash_top_tax_name'] == 'Drosophila melanogaster'
        assert str(row['sourmash_top_taxid']) == '11111'
        assert row['sourmash_order'] == 'Diptera'


def test_contig_sizes_give_every_contig_a_row():
    sizes = os.path.join(tempfile.mkdtemp(), 'asm.genome')
    with open(sizes, 'w') as f:
        f.write("contig_1\t1000\ncontig_2\t2000\ncontig_9\t500\n")
    with tempfile.TemporaryDirectory() as tmp:
        _, _, c_df = run(tmp, SAMPLE_SOURMASH_CSV, SAMPLE_ASSEMBLY_CSV, {'order': 'primates'}, contig_sizes=sizes)
        c = c_df.set_index('scaff')
        assert list(c_df['scaff']) == ['contig_1', 'contig_2', 'contig_9']
        assert c.loc['contig_9', 'sourmash_action'] == 'NA'
        assert c.loc['contig_9', 'sourmash_n_hits'] == 0
        assert c.loc['contig_9', 'sourmash_seq_len'] == 500
        assert c.loc['contig_2', 'sourmash_seq_len'] == 2000


def _fake_taxdump(d):
    """Tiny NCBI taxdump: 1 root > 2 Bacteria (blast name 'bacteria') > 976 Bacteroidota ('CFB group bacteria')
    > 1001 some bacterium; 7227 fly under 50557 Insecta ('insects'); 999 merged into 1001."""
    def w(name, rows):
        with open(os.path.join(d, name), 'w') as f:
            for r in rows:
                f.write("\t|\t".join(r) + "\t|\n")
    w('names.dmp', [("2", "Bacteria", "", "scientific name"), ("2", "bacteria", "", "blast name"),
                    ("976", "CFB group bacteria", "", "blast name"), ("50557", "insects", "", "blast name"),
                    ("7227", "Drosophila melanogaster", "", "scientific name")])
    w('taxidlineage.dmp', [("1001", "1 2 976 "), ("7227", "1 2759 50557 "), ("12345", "1 2759 ")])
    w('merged.dmp', [("999", "1001")])
    w('rankedlineage.dmp', [("7227", "Drosophila melanogaster", "", "Drosophila", "Drosophilidae", "Diptera",
                             "Insecta", "Arthropoda", "Metazoa", "Eukaryota")])
    mapping = os.path.join(d, 'map.tsv')
    with open(mapping, 'w') as f:
        f.write("# comment\n#BLAST-div\tGX-div\nbacteria\tprok:bacteria\nCFB group bacteria\tprok:CFB group bacteria\n"
                "insects\tanml:insects\n")
    return mapping


def test_gx_div_from_blast_names():
    d = tempfile.mkdtemp()
    mapping = _fake_taxdump(d)
    divs = stp.load_gx_divs(d, mapping, ['1001', '7227', '999', '12345', 'nan'])
    assert divs == {'1001': 'prok:CFB group bacteria', '7227': 'anml:insects',
                    '999': 'prok:CFB group bacteria', '12345': 'unkn:unknown'}


def test_gx_div_in_classification():
    d = tempfile.mkdtemp()
    mapping = _fake_taxdump(d)
    asm = ASM_HEADER + HOMO + DRO.replace('GCA_000003,7227', 'GCA_000003,7227')
    sm = "query_name,match_name,containment,jaccard,intersect_hashes\ncontig_1,GCA_000003,0.5,0.4,50\n"
    with tempfile.TemporaryDirectory() as tmp:
        _, _, c_df = run(tmp, sm, asm, {'order': 'primates'}, taxdump=d, div_mapping=mapping)
        assert c_df.set_index('scaff').loc['contig_1', 'sourmash_div'] == 'anml:insects'


def test_gx_div_missing_taxdump_files_is_empty():
    assert stp.load_gx_divs(tempfile.mkdtemp(), None, ['7227']) == {}


SPECIES_ASM = (ASM_HEADER.replace('taxid,', 'taxid,tax_name,')
               + "GCA_000001,9606,Homo sapiens,Homo sapiens,Homo,Hominidae,Primates,Mammalia,Chordata,Animalia\n"
               + "GCA_000004,63221,Homo sapiens neanderthalensis,Homo sapiens,Homo,Hominidae,Primates,Mammalia,Chordata,Animalia\n"
               + "GCA_000005,9598,Pan troglodytes,Pan troglodytes,Pan,Hominidae,Primates,Mammalia,Chordata,Animalia\n"
               + "GCA_000006,1052212,Homo sapiens gut symbiont,,,,,,,\n"
               + "GCA_000007,9606,Homo sapiens,,Homo,Hominidae,Primates,Mammalia,Chordata,Animalia\n")


def test_species_accessions_by_name_and_taxid():
    p = _tmp(SPECIES_ASM)
    try:
        df = stp.parse_assembly_database([p])
        # species column, or tax_name when species is empty (GCA_000007); never a symbiont named after the host
        assert stp.species_accessions(df, 'Homo sapiens') == {'GCA_000001', 'GCA_000004', 'GCA_000007'}
        # taxid resolved via the taxa table (no taxdump)
        assert stp.species_accessions(df, '9606') == {'GCA_000001', 'GCA_000004', 'GCA_000007'}
        assert stp.species_accessions(df, None) == set()
    finally:
        os.unlink(p)


def test_species_taxid_resolved_with_rankedlineage():
    d = tempfile.mkdtemp()
    _fake_taxdump(d)
    assert stp.species_name_for_taxid('7227', d) == 'Drosophila melanogaster'


def test_exclude_species_in_main():
    sm = ("query_name,match_name,containment,jaccard,intersect_hashes\n"
          "contig_1,GCA_000001,0.9,0.8,500\n"
          "contig_1,GCA_000005,0.5,0.4,300\n")
    with tempfile.TemporaryDirectory() as tmp:
        s_df, _, c_df = run(tmp, sm, SPECIES_ASM, {'order': 'primates'}, exclude_species='Homo sapiens')
        assert 'GCA_000001' not in s_df['assembly_accession'].values
        # the other Primates genome is still a target hit
        assert cls(c_df, 'contig_1') == 'target'
