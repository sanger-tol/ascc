#!/usr/bin/env bash
set -x
mkdir -p ascc_tiny_test_databases
cd ascc_tiny_test_databases
curl https://tolit.cog.sanger.ac.uk/test-data/resources/ascc/asccTinyTest_V2.tar.gz | tar xzf -
curl https://dp24.cog.sanger.ac.uk/ascc/diamond.dmnd -o diamond.dmnd
curl https://dp24.cog.sanger.ac.uk/blastn.tar.gz | tar xzf -
curl https://dp24.cog.sanger.ac.uk/ascc/accession2taxid.tar.gz | tar -xzf -
mkdir ncbi_taxdump
curl -L https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar -C ncbi_taxdump -xzf -
mkdir FCS_gx
wget -cq https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.taxa.tsv -O FCS_gx/all.taxa.tsv
wget -cq https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.gxi -O FCS_gx/all.gxi
wget -cq https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.gxs -O FCS_gx/all.gxs
wget -cq https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.meta.jsonl -O FCS_gx/all.meta.jsonl
wget -cq https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.blast_div.tsv.gz -O FCS_gx/all.blast_div.tsv.gz
mkdir busco_database
curl -L https://tolit.cog.sanger.ac.uk/test-data/resources/busco/blobtoolkit.GCA_922984935.2.2023-08-03.lineages.tar.gz | tar -C busco_database -xzf -
mkdir kraken2
curl -L https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/genome/db/kraken2.tar.gz | tar -C kraken2 -xzf -
mkdir NT_database
curl -L https://ftp.ncbi.nlm.nih.gov/blast/db/18S_fungal_sequences.tar.gz | tar -C NT_database -xzf -
mkdir pacbio_barcode
wget -O pacbio_barcode/SMRTbell_Barcoded_Adapter_Plate_3.0_bc2001-bc2096.fasta_.zip -c https://www.pacb.com/wp-content/uploads/SMRTbell_Barcoded_Adapter_Plate_3.0_bc2001-bc2096.fasta_.zip
cd pacbio_barcode
unzip SMRTbell_Barcoded_Adapter_Plate_3.0_bc2001-bc2096.fasta_.zip
mv SMRTbell_Barcoded_Adapter_Plate_3.0_bc2001-bc2096.fasta pacbio_adaptors.fa
rm -rf SMRTbell_Barcoded_Adapter_Plate_3.0_bc2001-bc2096.fasta_.zip __MACOSX
cd ../
mkdir diamond
wget -c https://dp24.cog.sanger.ac.uk/ascc/diamond.dmnd -O diamond/UP000000212_1234679_tax.dmnd
mkdir vecscreen
curl -L https://ftp.ncbi.nlm.nih.gov/blast/db/v4/16SMicrobial_v4.tar.gz | tar -C vecscreen -xzf -
