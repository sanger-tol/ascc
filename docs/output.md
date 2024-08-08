# sanger-tol/ascc: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [YamlInput](#yamlinput) -
- [Validate TaxID](#validate-taxid) -
- [Filter Fasta](#filter-fasta) -
- [GC Content](#gc-content) -
- [Generate Genome](#generate-genome) -
- [Trailing Ns Check](#trailing-ns-check) -
- [Get KMERS profile](#get-kmers-profile) -
- [Extract Tiara Hits](#extract-tiara-hits) -
- [Mito organellar blast](#mito-organellar-blast) -
- [Plastid organellar blast](#plastid-organellar-blast) -
- [Run FCS Adaptor](#run-fcs-adaptor) -
- [Run FCS-GX](#run-fcs-gx) -
- [Pacbio Barcode Check](#pacbio-barcode-check) -
- [Run Read Coverage](#run-read-coverage) -
- [Run Vecscreen](#run-vecscreen) -
- [Run NT Kraken](#run-nt-kraken) -
- [Nucleotide Diamond Blast](#nucleotide-diamond-blast) -
- [Uniprot Diamond Blast](#uniprot-diamond-blast) -
- [Create BTK dataset](#create-btk-dataset) -
- [Autofilter and check assembly](#autofilter-and-check-assembly) -
- [Generate samplesheet](#generate-samplesheet) -
- [Sanger-TOL BTK](#sanger-tol-btk) -
- [Merge BTK datasets](#merge-btk-datasets) -
- [ASCC Merge Tables](#ascc-merge-tables) -
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### YamlInput

<details markdown="1">
<summary>Output files</summary>

- `NA`

</details>

YamlInput parses the input yaml into channels for later use in the pipeline.


### Validate TaxID

<details markdown="1">
<summary>Output files</summary>

- `NA`

</details>

Validate TaxID scans through the taxdump to ensure that the input taxid is present in the nxbi taxdump.


### Filter Fasta

<details markdown="1">
<summary>Output files</summary>

- `filter/`
  `*filtered.fasta` - A fasta file that has been filtered for sequences below a given threshold.

</details>

By default scaffolds above 1.9Gb are removed from the assembly, as scaffolds of this size are unlikely to truely have contamination. There is also the issue that scaffolds larger than this use a significant amount of resources which hinders production environments.


### GC Content

<details markdown="1">
<summary>Output files</summary>

- `gc/`
  `*-GC_CONTENT.txt` - A text file describing the GC content of the input genome.

</details>

Calculating the GC content of the input genome.


### Generate Genome

<details markdown="1">
<summary>Output files</summary>

- `generate/`
  `*.genome` - An index-like file describing the input genome.

</details>

An index-like file containing the scaffold and scaffold length of the input genome.


### Trailing Ns Check

<details markdown="1">
<summary>Output files</summary>

- `trailingns/`
  `*_trim_Ns` - A text file containing a report of the Ns found in the genome.

</details>

A text file containing a report of the Ns found in the genome.


### Get KMERS profile

<details markdown="1">
<summary>Output files</summary>

- `get/`
  `*_KMER_COUNTS.csv` - A csv file containing kmers and their counts.

</details>

A csv file containing kmers and their counts.


### Extract Tiara Hits

<details markdown="1">
<summary>Output files</summary>

- `tiara/`
  `*.{txt,txt.gz}` - A text file containing classifications of potential contaminants.
  `log_*.{txt,txt.gz}` - A log of the tiara run.
  `*.{fasta,fasta.gz}` - An output fasta file.

</details>

Tiara ...


### Mito Organellar Blast

<details markdown="1">
<summary>Output files</summary>

- `blast/`
  `*.tsv` - A tsv file containing potential contaminants.

</details>

A BlastN based subworkflow used on the input genome to filter potential contaminants from the genome.


### Chloro Organellar Blast

<details markdown="1">
<summary>Output files</summary>

- `blast/`
  `*.tsv` - A tsv file containing potential contaminants.

</details>

A BlastN based subworkflow used on the input genome to filter potential contaminants from the genome.


### Run FCS Adaptor

<details markdown="1">
<summary>Output files</summary>

- `fcs/`
  `*.fcs_adaptor_report.txt` - A text file containing potential adaptor sequences and locations.
  `*.cleaned_sequences.fa.gz` - Cleaned fasta file.
  `*.fcs_adaptor.log` - Log of the fcs run.
  `*.pipeline_args.yaml` - Arguments to FCS Adaptor
  `*.skipped_trims.jsonl` - Skipped sequences

</details>

FCS Adaptor Identified potential locations of retained adaptor sequences from the sequencing run.


### Run FCS-GX

<details markdown="1">
<summary>Output files</summary>

- `fcs/`
  `*out/*.fcs_gx_report.txt` - A text file containing potential contaminant locations.
  `out/*.taxonomy.rpt` - Taxonomy report of the potential contaminants.

</details>

FCS-GX Identified potential locations of contaminant sequences.


### Pacbio Barcode Check

<details markdown="1">
<summary>Output files</summary>

- `filter/`
  `*_filtered.txt` - Text file of barcodes found in the genome.

</details>

Uses BlastN to identify where given barcode sequences may be in the genome.


### Run Read Coverage

<details markdown="1">
<summary>Output files</summary>

- `samtools/`
  `*.bam` - Aligned BAM file.
  `*_average_coverage.txt` - Text file containing the coverage information for the genome

</details>

Mapping the read data to the input genome and calculating the average coverage across it.


### Run Vecscreen

<details markdown="1">
<summary>Output files</summary>

- `summarise/`
  `*.vecscreen_contamination` - A text file containing potential vector contaminant locations.

</details>

Vecscreen identifies vector contamination in the input sequence.


### Run NT Kraken

<details markdown="1">
<summary>Output files</summary>

- `kraken2/`
  `*.classified{.,_}*'` - Fastq file containing classified sequence.
  `*.unclassified{.,_}*'` - Fastq file containing unclassified sequence.
  `*classifiedreads.txt` - A text file containing a report on reads which have been classified.
  `*report.txt` - Report of Kraken2 run.
- `get/`
  `*txt` - Text file containing lineage information of the reported meta genomic data.

</details>

Kraken assigns taxonomic labels to metagenomic DNA sequences and optionally outputs the fastq of these data.


### Nucleotide Diamond Blast

<details markdown="1">
<summary>Output files</summary>

- `diamond/`
  `*.txt` - A text file containing the genomic locations of hits and scores.
- `reformat/`
  `*text` - A Reformated text file continaing the full genomic location of hits and scores.
- `convert/`
  `*.hits` - A file containing all hits above the cutoff.

</details>

Diamond Blast is a sequence aligner for translated and protein sequences, here it is used do identify contamination usin the NCBI db


### Uniprot Diamond Blast

<details markdown="1">
<summary>Output files</summary>

- `diamond/`
  `*.txt` - A text file containing the genomic locations of hits and scores.
- `reformat/`
  `*text` - A Reformated text file continaing the full genomic location of hits and scores.
- `convert/`
  `*.hits` - A file containing all hits above the cutoff.

</details>

Diamond Blast is a sequence aligner for translated and protein sequences, here it is used do identify contamination usin the Uniprot db


### Create BTK dataset

<details markdown="1">
<summary>Output files</summary>

- `create/`
  `btk_datasets/` - A btk dataset folder containing data compatible with BTK viewer.
  `btk_summary_table_full.tsv` - A TSV file summarising the dataset.

</details>

Create BTK, creates a BTK_dataset folder compatible with BTK viewer.


### Autofilter and check assembly

<details markdown="1">
<summary>Output files</summary>

- `autofilter/`
  `autofiltered.fasta` - The decontaminated input genome.
  `ABNORMAL_CHECK.csv` - Combined FCS and Tiara summary of contamination.
  `assembly_filtering_removed_sequences.txt` - Sequences deemed contamination and removed from the above assembly.
  `fcs-gx_alarm_indicator_file.txt` - Contains text to control the running of Blobtoolkit.

</details>

Autofilter and check assembly returns a decontaminated genome file as well as summaries of the contamination found.


### Generate samplesheet

<details markdown="1">
<summary>Output files</summary>

- `generate/`
  `*.csv` - A CSV file containing data locations, for use in Blobtoolkit.

</details>

This produces a CSV containing information on the read data for use in BlobToolKit.


### Sanger-TOL BTK

<details markdown="1">
<summary>Output files</summary>

- `sanger/`
  `*_btk_out/blobtoolkit/${meta.id}*/` - The BTK dataset folder generated by BTK.
  `*_btk_out/blobtoolkit/plots/` - The plots for display in BTK Viewer.
  `*_btk_out/blobtoolkit/${meta.id}*/summary.json.gz` - The Summary.json file...
  `*_btk_out/busco/*` - The BUSCO results returned by BTK.
  `*_btk_out/multiqc/*` - The MultiQC results returned by BTK.
  `blobtoolkit_pipeline_info` - The pipeline_info folder.

</details>

Sanger-Tol/BlobToolKit is a Nextflow re-implementation of the snakemake based BlobToolKit pipeline and produces interactive plots used to identify true contamination and seperate sequence from the main assembly.


### Merge BTK datasets

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

- `ascc/`
  `*_contamination_check_merged_table.csv` - ....
  `*_contamination_check_merged_table_extended.csv` - ....
  `*_phylum_counts_and_coverage.csv` - A CSV report containing information on the hits per phylum and the coverage of the hits..

</details>

Merge Tables merged the summary reports from a number of modules inorder to create a single set of reports.


### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
