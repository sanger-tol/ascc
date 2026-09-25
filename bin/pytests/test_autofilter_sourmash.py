import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import autofilter  # noqa: E402

CLASSIFICATION = (
    "scaff,sourmash_top_tax_name,sourmash_top_taxid,sourmash_div,sourmash_action,sourmash_species\n"
    'ctg_1,"Wolbachia endosymbiont, strain X",1,prok:a-proteobacteria,EXCLUDE,x\n'
    "ctg_2,Hokovirus HKV1,2,virs:viruses,REVIEW,y\n"
    "ctg_3,Lagria hirta,3,anml:insects,NA,z\n"
)
NON_TARGET = "header,assembly_accession,taxa\nctg_1,GCF_1,1\nctg_4,GCF_2,2\n"


def _tmp(content):
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(content)
        return f.name


def test_classification_format_actions():
    p = _tmp(CLASSIFICATION)
    try:
        assert autofilter.get_sourmash_action_dict(p, "remove") == {"ctg_1": "EXCLUDE", "ctg_2": "REVIEW"}
        assert autofilter.get_sourmash_action_dict(p, "warn") == {"ctg_1": "WARN", "ctg_2": "REVIEW"}
    finally:
        os.unlink(p)


def test_legacy_non_target_format():
    p = _tmp(NON_TARGET)
    try:
        assert autofilter.get_sourmash_action_dict(p, "remove") == {"ctg_1": "EXCLUDE", "ctg_4": "EXCLUDE"}
        assert autofilter.get_sourmash_action_dict(p, "warn") == {"ctg_1": "WARN", "ctg_4": "WARN"}
    finally:
        os.unlink(p)


def test_empty_or_missing_file():
    p = _tmp("")
    try:
        assert autofilter.get_sourmash_action_dict(p) == {}
    finally:
        os.unlink(p)
    assert autofilter.get_sourmash_action_dict("/non/existent.csv") == {}
