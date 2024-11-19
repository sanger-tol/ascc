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

## NCBI nr Diamond database

Download the nr database protein FASTA files from the NCBI ftp server (`wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz`). The database building is similar to the Uniprot Diamond database building (as described at https://blobtoolkit.genomehubs.org/install/).
An example command for building the nr Diamond database looks like this:
```
diamond makedb --threads 16 --in ./nr/nr.gz -d 
./ncbi_taxonomy/proteins/nr --taxonmap ./ncbi_taxonomy/proteins/prot.accession2taxid.FULL --taxonnodes ./ncbi_taxonomy/proteins/taxdump/nodes.dmp --taxonnames ./ncbi_taxonomy/proteins/taxdump/names.dmp
```
The `prot.accession2taxid.FULL` file comes from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/.
The taxdump files come from ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz.


## NCBI accession2taxid

Download the files from the NCBI FTP server and uncompress them:

```
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/dead_wgs.accession2taxid.gz.md5 \
&& wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/dead_nucl.accession2taxid.gz.md5 \
&& wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz.md5 \
&& wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz.md5 \
&& wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/dead_nucl.accession2taxid.gz \
&& wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/dead_wgs.accession2taxid.gz \
&& wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz \
&& wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz \
&& gunzip *accession2taxid.gz
```

## FCS-GX database

Downloading instructions are at https://github.com/ncbi/fcs/wiki/FCS-GX-quickstart#download-the-fcs-gx-database.

## FCS-adaptor

The FCS-adaptor database is included in the FCS-adaptor installation, so it doesn't need to be downloaded separately.

## VecScreen

A FASTA file with the sequences for making a VecScreen database is included in the ASCC repository. It is the `vecscreen_adaptors_for_screening_euks.fa` file in the `assets` directory of this pipeline ([vecscreen_adaptors_for_screening_euks.fa](../assets/vecscreen_adaptors_for_screening_euks.fa)).

## PacBio barcodes

A FASTA file with the sequences of PacBio multiplexing barcodes is included in the ASCC repository. It is the `pacbio_adaptors.fa` file in the `assets` directory of this pipeline ([pacbio_adaptors.fa](../assets/pacbio_adaptors.fa)).
