[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/sanger-tol/ascc)

---

## Introduction

**sanger-tol/ascc** is a bioinformatics pipeline that is meant for detecting cobionts and contaminants in genome assemblies. ASCC stands for Assembly Screen for Cobionts and Contaminants. The pipeline aggregates tools such as BLAST, GC and coverage calculation, FCS-adaptor, FCS-GX, VecScreen, BlobToolKit, the BlobToolKit pipeline, Tiara, Kraken, Diamond BLASTX, and kmer counting and with kcounter+scipy. The main outputs are:

- A CSV table with taxonomic classifications of the sequences from the consitutent tools.
- A BlobToolKit dataset that can contain variables that are not present in BlobToolKit datasets produced by the BlobToolKit pipeline (https://github.com/sanger-tol/blobtoolkit) on its own. For example, ASCC can incorporate FCS-GX results into a BlobToolKit dataset.
- Individual report files for adapter, PacBio barcode and organellar contaminants.
  The only required input file for ASCC is the assembly FASTA file. Optional inputs are sequencing reads and organellar FASTA files. All individual components of the pipeline are optional, so it is possible to do lightweight runs with assemblies that have a simple composition of species and comprehensive runs with assemblies with complex composition.

![sanger-tol/ascc overview diagram](docs/images/ascc_overview_diagram.png)

1. Run a selection of processes from the list below (pick any that you think will be useful).

- FCS-GX
- FCS-adaptor
- VecScreen
- Tiara
- BlobToolKit Pipeline
- nt BLAST
- nr and Uniprot Diamond BLASTX
- GC and coverage calculation
- PacBio barcodes screen
- Organellar BLAST
- nt Kraken2
- kmer counting + dimensionality reduction

2. Postprocess the results of the previous step to produce summary files. What processes were run in the previous step determines what summary files can be generated. The possible outputs are:

- CSV table of sequence classification results
- BlobToolKit dataset
- CSV table of average coverage per phylum
- Adapter and organellar contamination report files

## Usage

> **Note**
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
> to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
> with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

-->

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run sanger-tol/ascc \
   -profile <docker/singularity/.../institute> \
   --input YAML \
   --outdir <OUTDIR> -entry SANGERTOL_ASCC --include ALL
```

> **Warning:**
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
> provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Credits

sanger-tol/ascc was written by [Eerik Aunin](https://github.com/eeaunin), [Damon Lee Pointon](https://github.com/DLBPointon), [James Torrance](https://github.com/jt8-sanger), [Ying Sims](https://github.com/yumisims) and [Will Eagles](https://github.com/weaglesBio). Pipeline development was supervised by [Shane A. McCarthy](https://github.com/mcshane) and [Matthieu Muffato](https://github.com/muffato).

We thank [Michael Paulini](https://github.com/epaule), Camilla Santos, [Noah Gettle](https://github.com/gettl008) and [Ksenia Krasheninnikova](https://github.com/ksenia-krasheninnikova) for testing the pipeline.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  sanger-tol/ascc for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
