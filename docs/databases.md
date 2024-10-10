# Databases used by the ASCC pipeline
This document describes the databases used by the pipeline.

## nt Kraken
For building the database, install Kraken2 if it is not already installed.
It can be installed using [conda](https://anaconda.org/bioconda/kraken2) with this command:
`conda install bioconda::kraken2`<br>
Commands for building the Kraken database:<br>
```
kraken2-build --download-library nt --db nt --threads 16<br>
kraken2-build --download-taxonomy --db /nt
```
Warning: building this database requires hundreds of gigabytes of memory.

## BUSCO5 database
Download from from:
https://busco-data.ezlab.org/v5/data/.

## NCBI taxdump, NCBI nt database and Uniprot Diamond database
Download and set up according to the instructions at https://blobtoolkit.genomehubs.org/install/.

## NCBI nr Diamond database
Download the nr database protein FASTA files from the NCBI ftp server (`wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz`) and build the database similarly to the Uniprot Diamond database, following the instructions at https://blobtoolkit.genomehubs.org/install/.

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