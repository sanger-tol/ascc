# sanger-tol/ascc: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## YAML input

### Full YAML

```yaml
scientific_name: scientific name of the assembled organism
taxid: NCBI taxonomy ID of the assembled species (or genus). Should be a numerical value, e.g. 352914. You can look up the TaxID for your species at https://ncbi.nlm.nih.gov/taxonomy
reads_path: path to a directory that contains gzipped reads
reads_type: determines which minimap2 preset will be used for read mapping. While minimap2 supports various read types (Illumina paired-end, PacBio CLR, PacBio HiFi, Oxford Nanopore), currently only "hifi" is implemented in this pipeline
pacbio_barcode_file: full path to the PacBio multiplexing barcode sequences database file. A FASTA file with known PacBio multiplexing barcode sequences is bundled with this pipeline, at "/ascc/assets/pacbio_adaptors.fa")
pacbio_barcode_names: comma separated list of names of PacBio multiplexing barcodes that were used in the sequencing of this sample. For example: "bc2008,bc2009". The barcode names exist in the barcode sequences database file ("/ascc/assets/pacbio_adaptors.fa")
kmer_length: kmer length for kmer counting (which is done using kcounter). Default: 7
dimensionality_reduction_methods: a comma separated list of methods for the dimensionality reduction of kmer counts. The available methods are the following: ["pca","umap","t-sne","isomap","lle_standard","lle_hessian","lle_modified","mds","se","random_trees","kernel_pca","pca_svd","autoencoder_sigmoid","autoencoder_linear","autoencoder_selu","autoencoder_relu","nmf"]. The default method is "pca". This field should be formatted as a YAML list, e.g. ["pca","random_trees"]
nt_database_path: path to the directory that contains the NCBI nt BLAST database. The database should have built-in taxonomy. Should end with a trailing slash
nt_database_prefix: prefix for the NCBI nt database. Default: "nt"
nt_kraken_database_path: path + prefix to the Kraken database made from NCBI nt database sequences
ncbi_accession_ids_folder: path to the directory with NCBI accession2taxid files (e.g. "/accession2taxid/"). Should end with a trailing slash
ncbi_taxonomy_path: path to NCBI taxdump directory (e.g. "/taxdump/"). Should end with a trailing slash
ncbi_ranked_lineage_path: path to NCBI ranked lineage file (e.g. "/taxdump/rankedlineage.dmp")
busco_lineages_folder: path to BUSCO 5 lineages directory. Should end with a trailing slash
busco_lineages: a comma separated list of BUSCO lineages that will be used in the sanger-tol/blobtoolkit pipeline run. For example: "diptera_odb10,insecta_odb10". Available lineages can be found at https://busco-data.ezlab.org/v5/data/lineages/
fcs_gx_database_path: path to the directory containing the FCS-GX database. Should end with a trailing slash
vecscreen_database_path: path to the FASTA file with adapter sequences for VecScreen ("/ascc/assets/vecscreen_adaptors_for_screening_euks.fa")
diamond_uniprot_database_path: path to a Diamond database made from Uniprot protein sequences ("uniprot_reference_proteomes_with_taxonnames.dmnd"). The database needs to have built-in taxonomy
diamond_nr_database_path: path to a Diamond database made from NCBI nr protein sequences ("nr.dmnd"). The database needs to have built-in taxonomy
seqkit_sliding: sliding window step size in bp, when sampling sequences for ASCC's built-in BLAST and Diamond processes. Default: 100000
seqkit_window: length of each sampled sequence in bp, when sampling sequences for ASCC's built-in BLAST and Diamond processes. Default: 6000
n_neighbours: n_neighbours setting for the kmers dimensionality reduction. This applies to the dimensionality reduction methods that have a n_neighbours parameter, such as UMAP. Default: 13
btk_yaml: path to a dummy YAML file that is provided with this pipeline, at "/ascc/assets/btk_draft.yaml". This is default and only serves to bypass GCA requirements of sanger-tol/blobtoolkit
```

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
Usage:
nextflow run sanger-tol/ascc \
    --input {INPUT YAML} \
    --outdir {OUTDIR} \
    [--include {COMMA SEPARATED LIST OF STEPS TO RUN}] \
    [--exclude {COMMA SEPARATED LIST OF STEPS TO EXCLUDE}] \
    [--organellar_include {COMMA SEPARATED LIST OF STEPS TO RUN}] \
    [--organellar_exclude {COMMA SEPARATED LIST OF STEPS TO EXCLUDE}] \
    -profile singularity
```

This will launch the pipeline with the `singularity` configuration profile. See below for more information about profiles.

Pipeline component options:

`--include`: comma-separated list of pipeline components to run on chromosomal DNA sequences (primary and haplotigs).<br>
`--exclude`: comma-separated list of pipeline components to exclude from running on chromosomal DNA sequences.<br>
`--organellar_include`: comma-separated list of pipeline components to run on organellar DNA sequences (mitochondrial and plastid).<br>
`--organellar_exclude`: comma-separated list of pipeline components to exclude from running on organellar DNA sequences.

Available pipeline components:
- `kmers`              : K-mer counting and dimensionality reduction analysis using kcounter, scikit-learn, and TensorFlow
- `tiara`              : Deep learning-based classification of sequences into prokaryotic and eukaryotic origin using Tiara
- `coverage`           : Analysis of sequence coverage using minimap2-based read mapping
- `nt_blast`          : Nucleotide BLAST search against NCBI nt database for taxonomic classification
- `nr_diamond`        : DIAMOND BLASTX search against NCBI non-redundant protein database
- `uniprot_diamond`   : DIAMOND BLASTX search against UniProt database
- `kraken`            : Taxonomic classification using Kraken2 against NCBI nt database
- `fcs-gx`            : NCBI's FCS-GX (Foreign Contamination Screen - Genome Cross-Check) for contamination detection
- `fcs-adaptor`       : NCBI's FCS-Adaptor (Foreign Contamination Screen - Adaptor Screen) for adapter contamination
- `vecscreen`         : NCBI's vector contamination screening
- `btk_busco`         : BlobToolKit Pipeline (sequence classification using BUSCO, Diamond and BLAST)
- `pacbio_barcodes`   : Detection of PacBio barcode contamination
- `organellar_blast`  : BLAST-based detection of organellar sequences
- `autofilter_assembly`: Automated assembly filtering (requires `tiara` and `fcs-gx`)
- `ALL`               : Run all available components
- `NONE`              : Run no components

Dependencies:
- `autofilter_assembly` requires both `tiara` and `fcs-gx` to be run first

Outputs:
- Results are collected as BlobToolKit datasets and CSV tables
- Adapter and organellar contamination reports are provided as text files


### Example usage

#### Basic run with essential components
```
nextflow run sanger-tol/ascc --input config.yaml --outdir results --include tiara,coverage,nt_blast --organellar_include nt_blast,coverage -profile singularity
```

#### Comprehensive analysis
```
nextflow run sanger-tol/ascc --input config.yaml --outdir results --include kmers,tiara,coverage,nt_blast,nr_diamond,kraken,fcs-gx,btk_busco --organellar_include nt_blast,coverage -profile singularity
```

#### Run everything except specific components
```
nextflow run sanger-tol/ascc --input config.yaml --outdir results --include ALL --exclude vecscreen,pacbio_barcodes --organellar_include ALL -profile singularity
```

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull sanger-tol/ascc
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [sanger-tol/ascc releases page](https://github.com/sanger-tol/ascc/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> 💡 If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
