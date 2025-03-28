# sanger-tol/ascc: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

### Processes that produce the main outputs:

- [Trailing Ns Check](#trailing-ns-check)
- [Mito Organellar Blast](#mito-organellar-blast) -
- [Plastid organellar blast](#plastid-organellar-blast) -
- [Run FCS Adaptor](#run-fcs-adaptor) -
- [Pacbio Barcode Check](#pacbio-barcode-check) -
- [Run Read Coverage](#run-read-coverage) -
- [Run VecScreen](#run-vecscreen) -
- [Create BTK Dataset](#create-btk-dataset) -
- [Autofilter and Check Assembly](#autofilter-and-check-assembly) -
- [Sanger-TOL BTK](#sanger-tol-btk) -
- [Merge BTK datasets](#merge-btk-datasets) -
- [ASCC Merge Tables](#ascc-merge-tables) -
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Processes that produce intermediate outputs:

- [YamlInput](#yamlinput) -
- [Generate samplesheet](#generate-samplesheet) -
- [Validate TaxID](#validate-taxid) -
- [Generate Genome](#generate-genome) -
- [Filter Fasta](#filter-fasta) -
- [GC Content](#gc-content) -
- [Get kmers profile](#get-kmers-profile) -
- [Extract Tiara Hits](#extract-tiara-hits) -
- [Run FCS-GX](#run-fcs-gx) -
- [Run nt Kraken](#run-nt-kraken) -
- [nr Diamond BLASTX](#nr-diamond-blastx) -
- [Uniprot Diamond BLASTX](#uniprot-diamond-blastx) -

## Main outputs

### Trailing Ns Check

<details markdown="1">
<summary>Output files</summary>

- `trailingns/`
  `*_trim_Ns` - A text file containing a report of trailing Ns found in the genome.

</details>

A text file containing a report of trailing Ns found in the genome. Trailing Ns are when a nucleotide sequence starts or ends with Ns instead of A, G, C or T nucleotides. It is advisable to trim off the trailing Ns from sequences in the assembly. If the sequence remaining after trimming is shorter than 200 bp, the script recommends removing it from the assembly.

### Mito Organellar Blast

<details markdown="1">
<summary>Output files</summary>
- `organelle/`
  `*-mitochondrial_genome.contamination_recommendation` - A file that contains the names of sequences that are suspected mitochondrial contaminants in the nuclear DNA assembly, tagged as either "REMOVE" or "Investigate" depending on the BLAST hit alignment length and percentage identity. The file is empty if there are no suspected mitochondrial contaminants.
</details>

This subworkflow uses BLAST against a user-provided mitochondrial sequence to detect leftover organellar sequences in the assembly file that should contain only chromosomal DNA sequences. A BLAST nucleotide database is made from the user-provided organellar sequence. BLAST with the chromosomal DNA assembly is then ran against this database with the following settings: -task megablast -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 0.0001 -perc_identity 80 -soft_masking true. The BLAST results are filtered to keep only hits with alignment length that is at least 200 bp.
Depending on the alignment length and percentage identity, the script can recommend an action for dealing with the putative organellar sequence: either "REMOVE" or "Investigate".

### Plastid Organellar Blast

<details markdown="1">
- `organelle/`
  `*-plastid_genome.contamination_recommendation` - A file that contains the names of sequences that are suspected plastid contaminants in the nuclear DNA assembly, tagged as either "REMOVE" or "Investigate" depending on the BLAST hit alignment length and percentage identity. The file is empty if there are no suspected mitochondrial contaminants.

</details>

This subworkflow uses BLAST against a user-provided plastid sequence to detect leftover organellar sequences in the assembly file that should contain only chromosomal sequences. The method is the same as in the Mito Organellar Blast part.

### Run FCS-adaptor

<details markdown="1">
<summary>Output files</summary>

- `fcs/`
  `*.fcs_adaptor_report.txt` - A text file containing potential adaptor sequences and locations.
  `*.cleaned_sequences.fa.gz` - Cleaned FASTA file.
  `*.fcs_adaptor.log` - Log of the FCS-adaptor run.
  `*.pipeline_args.yaml` - Arguments to FCS-adaptor
  `*.skipped_trims.jsonl` - Skipped sequences

</details>

FCS-adaptor (https://github.com/ncbi/fcs) is NCBI software for detecting adapter contamination in genome assemblies. FCS-adaptor uses a built-in database of adapter sequences, provided by NCBI. The FCS-adaptor report shows identified potential locations of retained adapter sequences from the sequencing run.

### Pacbio Barcode Check

<details markdown="1">
<summary>Output files</summary>

- `filter/`
  `*_filtered.txt` - Text file log of PacBio barcode sequences found in the genome. The file is empty if no contamination was found.

</details>

Uses BlastN to identify retained PacBio multiplexing barcode contamination in the assembly. The PacBio multiplexing barcode sequences are stored as the pacbio_adaptors.fa file in the assets directory of this pipeline.

### Run Read Coverage

<details markdown="1">
<summary>Output files</summary>
  `*.bam` - BAM file with aligned reads.
  `*_average_coverage.txt` - Text file containing the coverage information for the genome

</details>

Mapping the read data to the input genome with minimap2 (https://github.com/lh3/minimap2) and calculating the average coverage per sequence. The reads used for mapping can be PacBio HiFi reads or paired end Illumina reads.

### Run VecScreen

<details markdown="1">
<summary>Output files</summary>

- `summarise/`
  `*.vecscreen_contamination` - A text file containing potential vector contaminant locations. The file is empty if no potential contaminants were found.

</details>

VecScreen (https://www.ncbi.nlm.nih.gov/tools/vecscreen/) is a tool for detecting adapter and vector contamination in genome assemblies. It is an older tool than FCS-adaptor. Its advantage over FCS-adaptor is that it can use a custom database of contaminant sequences made by the user, whereas FCS-adaptor comes with its built-in database.

### Create BTK Dataset

<Details markdown="1">
<summary>Output files</summary>

- `create/`
  `btk_datasets/` - A BlobToolKit (https://blobtoolkit.genomehubs.org) dataset folder containing data compatible with BTK viewer (https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/opening-a-dataset-in-the-viewer/).
  `btk_summary_table_full.tsv` - A TSV file summarising the contents of the BlobToolKit dataset. This file is created using the `blobtools filter --table` command of BlobToolKit.

</details>

Creates a BlobToolKit dataset folder compatible with BlobToolKit viewer (https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/opening-a-dataset-in-the-viewer/). The BlobToolKit dataset create by ASCC can contain much more variables than what the BlobToolKit pipeline (https://github.com/sanger-tol/blobtoolkit) produces.

### Autofilter and Check Assembly

<details markdown="1">
<summary>Output files</summary>

- `autofilter/`
  `autofiltered.fasta` - The decontaminated input genome. The decontamination is based on the results of FCS-GX.
  `ABNORMAL_CHECK.csv` - Combined FCS-GX and Tiara summary of contamination.
  `assembly_filtering_removed_sequences.txt` - Sequences deemed contamination by FCS-GX (labelled with the EXCLUDE tag by FCS-GX) and removed from the above assembly.
  `fcs-gx_alarm_indicator_file.txt` - Contains text to control the running of BlobToolKit pipeline. If enough contamination is found by FCS-GX, an alarm is triggered to switch on the running of BlobToolKit pipeline.

</details>

Autofilter and check assembly returns a decontaminated genome file as well as summaries of the contamination found.

### Sanger-TOL BTK

<details markdown="1">
<summary>Output files</summary>

- `sanger-tol-btk/`
  `*_btk_out/blobtoolkit/${meta.id}*/` - The BlobToolKit dataset folder generated by the sanger-tol/blobtoolkit pipeline.
  `*_btk_out/blobtoolkit/plots/` - BlobToolKit plots as PNG images, exported from the BlobToolKit dataset using blobtk (https://pypi.org/project/blobtk/).
  `*_btk_out/blobtoolkit/${meta.id}*/summary.json.gz` - The summary.json.gz file of the BlobToolKit dataset. It contains assembly metrics such as
  `*_btk_out/busco/*` - The BUSCO results returned by BlobToolKit.
  `*_btk_out/multiqc/*` - The MultiQC results returned by BlobToolKit.
  `blobtoolkit_pipeline_info` - The pipeline_info folder.

</details>

Sanger-Tol/BlobToolKit (https://github.com/sanger-tol/blobtoolkit) is a Nextflow re-implementation of the Snakemake based BlobToolKit pipeline (https://github.com/blobtoolkit/pipeline) and produces interactive plots used to identify contamination or cobionts and separate these sequences from the main assembly.

### Merge BTK Datasets

<details markdown="1">
<summary>Output files</summary>

- `merge/`
  `merged_datasets` - A BTK dataset.
  `merged_datasets/btk_busco_summary_table_full.tsv` - A TSV file containing a summary of the btk busco results.

</details>

This module merged the Create_btk_dataset folder with the Sanger-tol BTK dataset to create one unified dataset for use with btk viewer.

### ASCC Merge Tables

<details markdown="1">
<summary>Output files</summary>

- `ASCC-main-output/`
`*_contamination_check_merged_table.csv` - A CSV table that contains the results of most parts of the pipeline (GC content, coverage, Tiara, Kraken, kmers dimensionality reduction, Diamond, BLAST, FCS-GX, BlobToolKit pipeline) for each sequence in the input assembly file.
If a set of prerequisite steps have been run (nt BLAST, nr Diamond, Uniprot Diamond, read mapping for coverage calculation, Tiara, nt Kraken and the creation of a BlobToolKit dataset), the pipeline tries to put together a phylum level combined classification of the input sequences. It first uses BlobToolKit's `bestsum_phylum`, then fills the gaps (caused by `no-hit` sequences) with results from Tiara and then the remaining gaps are filled with results from nt Kraken. The combined classification is in the `merged_classif` column. The `merged_classif_source` column says which tool's output the classification for each sequence is based on. The automated classification usually has some flaws in it but is still useful as a starting point for determining the phyla that the input sequences belong to.
 `*_phylum_counts_and_coverage.csv` - A CSV report containing information on the hits per phylum and the average coverage per phylum. This file can only be generated if the`merged_classif` variable has been produced in the `*_contamination_check_merged_table.csv` table, as described above.
</details>

Merge Tables merged the summary reports from a number of modules in order to create a single set of reports.

### Pipeline Information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

## Intermediate outputs

These files are produced by the pipeline's modules but they are stay in Nextflow's work directory and are not included on their own in the final output.

### YamlInput

<details markdown="1">
<summary>Output files</summary>

- `NA`

</details>

YamlInput parses the input YAML file into channels for later use in the pipeline.

### Validate TaxID

<details markdown="1">
<summary>Output files</summary>

- `NA`

</details>

Validate TaxID scans through the NCBI taxdump file to ensure that the taxonomy ID (taxID) provided by the user is present in the taxdump. The taxdump originates from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/. The taxID might be absent in the taxdump either because the user has provided a faulty value for taxID or the taxdump is out of date.

### Filter FASTA

<details markdown="1">
<summary>Output files</summary>

`*filtered.fasta` - A FASTA file that has been filtered to keep sequences below a given threshold of length.

</details>

By default scaffolds above 1.9 Gb are removed from the assembly, as scaffolds of this size are unlikely to truely have contamination. There is also the issue that scaffolds larger than this use a significant amount of resources which hinders production environments. Furthermore, [FCS-GX](https://github.com/ncbi/fcs) does not work with sequences larger than 2 Gb.

### GC Content

<details markdown="1">
<summary>Output files</summary>

`*-GC_CONTENT.txt` - A tab separated table describing the GC content of the input genome. The first column contains the sequence names and the second column contains the GC content of each sequence. The GC content is expressed as a fraction: number of G and C nucleotides in the sequence divided by the number of all nucleotides in the sequence.

</details>

Calculating the GC content of each sequence in the input genome.

### Generate Genome

<details markdown="1">
<summary>Output files</summary>

`*.genome` - An index-like file describing the input genome.

</details>

An index-like file containing the scaffold and scaffold length of the input genome.

### Get kmers Profile

<details markdown="1">
<summary>Output files</summary>

`*_KMER_COUNTS.csv` - A CSV file containing the counts of kmers (by default: 7mers) in each sequence in the assembly.
`KMERS_dim_reduction_embeddings_combined.csv` - A CSV file with the results of dimensionality reduction of kmer counts. The dimensionality reduction embeddings help to separate sequences in the assembly by their origin (sequences originating from the same species likely appear close together in an embedding). When setting up a run, the user can choose multiple methods for dimensionality reduction.

</details>

A CSV file containing the counts of kmers (by default: 7mers) in each sequence in the assembly. Also, a file with the results of dimensionality reduction of kmer counts.
The following dimensionality reduction methods are available: PCA (principal component analysis), kernel PCA, PCA with SVD (singular value decomposition) solver, UMAP (uniform manifold approximation and projection), t-SNE (t-distributed stochastic neighbor embedding), LLE (locally linear embedding), MDS (multidimensional scaling), SE (spectral embedding), random trees, autoencoder and NMF (non-negative matrix factorisation).
The first two dimensions of the dimensionality reduction embeddings are used as the x and y coordinate when visualising the results in BlobToolKit.

### Extract Tiara Hits

<details markdown="1">
<summary>Output files</summary>

`TIARA.txt` - A text file containing classifications of the input DNA sequences. Each sequence gets assigned one label out of these: archaea, bacteria, prokarya, eukarya, organelle and unknown).
`log_*.{txt}` - A log of the Tiara run.

</details>

Tiara (https://github.com/ibe-uw/tiara) uses a neural network to classify DNA sequences.

### Run FCS-GX

<details markdown="1">
<summary>Output files</summary>

- `fcs/`
  `*out/*.fcs_gx_report.txt` - A text file containing potential contaminant locations.
  `out/*.taxonomy.rpt` - Taxonomy report of the potential contaminants.

</details>

FCS-GX (https://github.com/ncbi/fcs) is NCBI software that detects contaminants in genome assemblies using a cross species aligner. It uses its own database, provided by NCBI.

### Run nt Kraken

<details markdown="1">
<summary>Output files</summary>

`*.kraken2.classifiedreads.txt` - A text file containing classifications for each input DNA sequence, generated by Kraken2.
`*.kraken2.report.txt` - Summary of the Kraken2 run, generated by Kraken2.
`_nt_kraken_lineage_file.txt` - Kraken2 lineages for each input DNA sequence, reformatted as a CSV table to make it possible to merge this information into a table that contains sequence classifications from other tools, e.g. BLAST and Diamond.

</details>

Kraken (https://github.com/DerrickWood/kraken2) assigns taxonomic labels to input DNA sequences based on comparing them to a database of kmers. ASCC uses a Kraken database made from the sequences of the NCBI nt database. The FASTA sequences of NCBI nt database are available at https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/.

### nr Diamond BLASTX

<details markdown="1">
<summary>Output files</summary>

`*.txt` - A tabular text file containing the raw output of running Diamond BLASTX with sampled chunks of the assembly. The file contains BLASTX hits and scores/ Format: outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms sphylums salltitles
`full_coords.tsv`: A tabular text file containing the results from Diamond BLASTX where the coordinates of the BLASTX of chunks of assembly have been converted to coordinates in the full sequences of the assembly.
`*_diamond_blastx_top_hits.csv` - A file containing Diamond BLASTX top hits for each sequence in the input assembly file.
`*_diamond_outfmt6.tsv` - the `full_coords.tsv` file reformatted to make it compatible with BlobToolKit, so that the hits in it can be added to a BlobToolKit dataset. Format: outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

</details>

Diamond (https://github.com/bbuchfink/diamond) is a sequence aligner for protein sequences and translated nucleotide sequences. Here it is used to identify contamination using the NCBI nr database. The FASTA sequences of NCBI nr database are available at https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/.

### Uniprot Diamond BLASTX

<details markdown="1">
<summary>Output files</summary>

`*.txt` - A tabular text file containing the raw output of running Diamond BLASTX with sampled chunks of the assembly. The file contains BLASTX hits and scores/ Format: outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms sphylums salltitles
`full_coords.tsv`: A tabular text file containing the results from Diamond BLASTX where the coordinates of the BLASTX of chunks of assembly have been converted to coordinates in the full sequences of the assembly.
`*_diamond_blastx_top_hits.csv` - A file containing Diamond BLASTX top hits for each sequence in the input assembly file.
`*_diamond_outfmt6.tsv` - the `full_coords.tsv` file reformatted to make it compatible with BlobToolKit, so that the hits in it can be added to a BlobToolKit dataset. Format: outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

</details>
Diamond (https://github.com/bbuchfink/diamond) is a sequence aligner for protein sequences and translated nucleotide sequences. Here it is used to identify contamination using the Uniprot database.

### Generate Samplesheet

<details markdown="1">
<summary>Output files</summary>

- `generate/`
  `*.csv` - A CSV file containing data locations, for use in BlobToolkit.

</details>

This produces a CSV containing information on the read data for use in BlobToolKit.
