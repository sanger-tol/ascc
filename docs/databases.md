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
