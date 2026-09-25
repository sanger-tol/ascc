# Databases used by the ASCC pipeline

This document describes the databases used by the pipeline.

## nt Kraken

For building the database, install Kraken2 if it is not already installed.
It can be installed using [conda](https://anaconda.org/bioconda/kraken2) with this command:

```
conda install bioconda::kraken2
```

Commands for building the Kraken database:<br>

```
kraken2-build --threads 16 --download-taxonomy --db nt
kraken2-build --threads 16 --download-library nt --db nt
kraken2-build --build --threads 16 --db nt
```

Warning: building this database requires hundreds of gigabytes of memory.
See their Manual at: https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown

## BUSCO5 database

Download from from:
https://busco-data.ezlab.org/v5/data/.

## NCBI taxdump, NCBI nt database and Uniprot Diamond database

Download and set up according to the instructions at https://blobtoolkit.genomehubs.org/install/.

**Important Note**: The NCBI nt BLAST database must have built-in taxonomy for the pipeline to work correctly. When building the BLAST database, make sure to include taxonomy information using the `-parse_seqids` and `-taxid_map` options with `makeblastdb`. Alternatively, you can download the pre-built BLAST database from NCBI which already includes taxonomy information.

## NCBI nr Diamond database

Download the nr database protein FASTA files from the NCBI ftp server (`wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz`) and build the database similarly to the Uniprot Diamond database, following the instructions at https://blobtoolkit.genomehubs.org/install/.

## FCS-GX database

Downloading instructions are at https://github.com/ncbi/fcs/wiki/FCS-GX-quickstart#download-the-fcs-gx-database.

## Sourmash databases

Sourmash searches (`sourmash scripts multisearch`, branchwater plugin) run against the databases listed in the CSV passed with `--sourmash_db_config` (production: [sourmash_databases_production.csv](../assets/sourmash_databases_production.csv)). Columns: `name,path,k_available,k_for_search,s,assembly_taxa_db`, where `path` is a RocksDB index and `assembly_taxa_db` is a CSV mapping each assembly accession to its NCBI lineage.

At Sanger the databases are in `/lustre/scratch124/tol/teams/grit/users/dz11/tax_classification_db/sourmash_db/`:

| Name                 | Genomes | k   | scaled | Source                                                     |
| -------------------- | ------- | --- | ------ | ---------------------------------------------------------- |
| vertebrate           | 1,540   | 51  | 200    | NCBI eukaryote "reference genome" assemblies (August 2025) |
| invertebrate         | 3,081   | 51  | 200    | as above                                                   |
| plant_fungi_protazoa | 2,560   | 51  | 200    | as above                                                   |
| gtdb_v226            | 143,614 | 31  | 200    | GTDB r226 species representatives                          |
| ncbi_viruses         | 230,973 | 21  | 50     | prebuilt sourmash NCBI viruses database (January 2025)     |

The references were not screened for contamination when these databases were built. An inventory of every genome in them, with its lineage and its NCBI assembly status and flags (`excluded_from_refseq`), is in `/lustre/scratch124/tol/teams/grit/users/dz11/tax_classification_db/db_audit_2026-09/` (`sourmash_db_genomes.tsv`, built by `build_db_genomes_table.py` from the manifests of the signature archives the indexes were built from).

[sourmash_exclude_accessions.tsv](../assets/sourmash_exclude_accessions.tsv) lists the assemblies whose hits should be ignored:

- NCBI flags them as contaminated, chimeric, mixed culture, misassembled or genome length too large;
- or they were found to contain foreign sequence during ASCC testing (the `reason` column says which).

These assemblies are still in the indexes; they are removed when the search results are parsed.

## FCS-adaptor

The FCS-adaptor database is included in the FCS-adaptor installation, so it doesn't need to be downloaded separately.

## VecScreen

A FASTA file with the sequences for making a VecScreen database is included in the ASCC repository. It is the `vecscreen_adaptors_for_screening_euks.fa` file in the `assets` directory of this pipeline ([vecscreen_adaptors_for_screening_euks.fa](../assets/vecscreen_adaptors_for_screening_euks.fa)).

VecScreen requires a BLAST V4 database as input, we can generate this with the above file use the following.

```
makeblastdb -in vecscreen_adaptors_for_screening_euks.fa -parse_seqids -blastdb_version 4 -dbtype nucl
```

To use this database, point the `vecscreen_database_path` variable in the input YAML file of the pipeline run to the directory that contains this BLAST database. Use the name of the directory for `vecscreen_database_path`, without using the name of the database files. E.g. `/path/to/my/database/files/vecscreen_database/`.

## PacBio barcodes

A FASTA file with the sequences of PacBio multiplexing barcodes is included in the ASCC repository. It is the `pacbio_adaptors.fa` file in the `assets` directory of this pipeline ([pacbio_adaptors.fa](../assets/pacbio_adaptors.fa)).
