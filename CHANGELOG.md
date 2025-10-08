# sanger-tol/ascc: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.6.0] - Red Notebook [##/##/2025]

THIS IS STILL AN IN-DEVELOPMENT PROJECT SO THERE MAY BE BUGS.

Release 10 of sanger-tol/ascc, addition of a report generator.

### `Added`
- Reporting for the pipeline output into a more human readable output (html).
- Updated the naming of various outputs inorder to standardise them.

### `Dependencies`

| Module                        | Old Version                           | New Versions                          |
| ----------------------------- | ------------------------------------- | ------------------------------------- |
| AUTOFILTER_AND_CHECK_ASSEMBLY | abnormal_contamination_check.py:1.1.0 | abnormal_contamination_check.py:1.2.0 |

## [0.5.2] - Red Spider-Boat (H2) [06/10/2025]

THIS IS STILL AN IN-DEVELOPMENT PROJECT SO THERE MAY BE BUGS.

Release 9 of sanger-tol/ascc, a modification to FCS_ADAPTOR and configs.

### `Added`

- Patch to `FCS_ADAPTOR` to avoid the use of `/tmp`
- Update `SANGER_TOL_BLOBTOOLKIT` in `base.config` to use 1200.MB rather than `process_low`'s 12.GB
- Updated the script `abnormal_contamination_check.py` in the module `AUTOFILTER_AND_CHECK_ASSEMBLY`
  - This adds reporting for the number of `REVIEW/INFO` tags output by FCS as an alarm paramter to trigger SANGER_TOL_BLOBTOOLKIT
- Updated ro-crate, tests, and CHANGELOG.
- Now the pipeline is quite stable in production, there is the aim to once again start collecting resource statistics.
  - Updated `trace` scope output for ASCC and BLOBTOOLKIT (via `assets/btk_config_files/btk_trace.config`)
- Updating `SANGER_TOL_BTK` to 0.9.0 (Scyther)

### `Dependencies`

| Module                        | Old Version                           | New Versions                          |
| ----------------------------- | ------------------------------------- | ------------------------------------- |
| AUTOFILTER_AND_CHECK_ASSEMBLY | abnormal_contamination_check.py:1.1.0 | abnormal_contamination_check.py:1.2.0 |
| SANGER_TOL_BTK                | 0.8.0                                 | 0.9.0                                 |

## [0.5.1] - Red Spider-Boat (H1) [29/10/2025]

THIS IS STILL AN IN-DEVELOPMENT PROJECT SO THERE MAY BE BUGS.

Release 8 of sanger-tol/ascc, focusing on a module update and some minor structural updates.

### `Added`

- `SAMTOOLS_FAIDX` now outputs a `.sizes` file so `CUSTOM_CHROMSIZES` has been removed.
- Remove Legacy `GrabFiles` function (#163).
- Updating modules to most recent version available on nf-core (#163, #162).
- Updated Snapshots to reflect version changes (#163)

### `Fixed`

- Removed nf-core modules which arn't actually in use (#164)
- Updated version outputs from modules previously reporting `null` (#164)

### `Bugs`

- Currently, Blobtoolkit will _not_ run if there is no autofilter output channel.
- tiara, fcsgx, autofilter must always be activated.

### `Dependencies`

| Module               | Old Version               | New Versions |
| -------------------- | ------------------------- | ------------ |
| SAMTOOLS_FAIDX       | 1.21.1                    | 1.22.1       |
| SAMTOOLS_SORT        | 1.21.1                    | 1.22.1       |
| SAMTOOLS_MERGE       | 1.21.1                    | 1.22.1       |
| SAMTOOLS_DICT        | 1.21.1                    | 1.22.1       |
| SAMTOOLS_DEPTH       | 1.21.1                    | 1.22.1       |
| KRAKEN2_KRAKEN2      | 2.1.4                     | 2.1.6        |
| BLASTN               | 2.15.0                    | 2.16.0       |
| BLAST_MAKEBLASTDB    | 2.15.0                    | 2.16.9       |
| DIAMOND_BLASTX       | 2.1.8                     | 2.1.12       |
| GNU_SORT             | 9.3                       | 9.5          |
| SEQKIT_SLIDER        | 2.8.1                     | 2.8.0        |
| CUSTOM_GETCHROMSIZES | samtools:1.21--h50ea8bc_0 | REMOVED      |

## [0.5.0] - Red Spider-Boat [05/10/2025]

THIS IS STILL AN IN-DEVELOPMENT PROJECT SO THERE MAY BE BUGS.

Release 7 of sanger-tol/ascc, focusing on the 3.3.2 template upgrade and stability for sanger production.

### `Notes`

- If running the pipeline in `--profile singularity` and you are crashing with a generic error, try using `export NXF_SINGULARITY_NEW_PID_NAMESPACE=false`.

### `Added`

- Template update to 3.3.2 (#155)
- Added pipeline-level nf-test testing which is now running as standard CI.
- Corrected versioning in the .nextflow.log.
- Minor updates to the base.config (#158)
- singularity pid setting is now false.
- Param to expose FILTER_FASTA ext.cutoff and enforce min/max values.
- The FCSGX module has been heavily patched when using `--profile production`, this is to support `module` and `modulecmd`.
- The `--production` profile now contains a FCSGX module override linked to the above.
  - Resources have been included otherwise FCSGX will use the nextflow defaults.
  - Queue for FCSGX has been changed to `oversubscribed` in the `production.config` (#162)
- Changes to the resource allocation to improve support for large genomes. Changes are for the modules:
  - BLAST_BLASTN
  - BLAST_BLASTN_MOD
  - DIAMOND_BLASTX
  - TIARA_TIARA
  - MINIMAP2_ALIGN_SE
- Above changes also seperate out the logic for organellar genomes as they do not need the more complicated resource additions.
- Removed the BTK config, this has started to cause a number of crashes to the nested blobtoolkit pipeline.
- Remove unnecessary SAMTOOLS_INDEX from RUN_READ_COVERAGE.

### `Fixed`

- Spelling mistakes... again.
- singularity pid env change is down to an issue with Singularity and FCSGX, multiple instances of the tool accessing the same DB files causes crashes.
- Correct value of 100Mb to 1Gb as the ext.cutoff for FILTER_FASTA
- GENERATE_SAMPLESHEET was only taking into account 1 read file for BTK rather than all provided read files.
- GENERATE_SAMPLESHEET was not using the `reads_type` variable.
- Changed instances of `projectDir` to `launchDir` for safety (In tests).
- Map pattern in ESSENTIAL_JOBS has been updated to reduce re-writing all values.
- Corrected references to `withName: '*:PACBIO_BARCODE_CHECK:BLAST_BLASTN'`.
- Updated resource notation to use e notation (#134)

### `Bugs`

- Currently, Blobtoolkit will _not_ run if there is no autofilter output channel.
- tiara, fcsgx, autofilter must always be activated.

## [0.4.0] - Red Mouse [ ##/06/2024]

THIS IS STILL AN IN-DEVELOPMENT PROJECT SO THERE MAY BE BUGS.

Release 6 of sanger-tol/ascc, focusing on template upgrade and stability for sanger production.

### `Added`

- Return free-disk-space to nf-test CI runners (suggestion by @prototaxites)
- Updated test files to post-datacentre crash (Only affects internal sanger users).
- Test and Production (this is sanger specific, please change if you want to use the same style) configs have been updated. #106 #111
- Update to Organellar Blast subworkflow to include organellar name in output.
- Updated some scripts for logging and linting.
- The AUTOFILTER_AND_CHECK_ASSEMBLY module has been added to the ORGANELLAR subworkflow, done to ensure that sanger-tol/ascc matches previously reported statistics.
- FCSGX Module now includes a timed backoff for when it fails, this is aimed at stopping random crashes.
- Re-organised some of the files in the assets folder and updated configs to reflect.
- Added a `btk_pipeline.config` file in `assets/btk_config_files` to modify btk process resource requirements. Currently contains an alternate requirement for `BLASTN_TAXON`.
- Added most of the output files to the emit of the major subworkflow (GENOMIC and ORGANELLAR) this is setup for future version nf-tools which will mandate output files are treated like this.
- Added `--fcs_override` and `--fcs_override_samplesheet` to allow the pipeline to accept externally run FCS-GX results. These results must be filteres and parsed as inside the pipeline. A wrapper script is provided as `bin/ascc_fcsgx_wrapper.py`.
- Support for `HAP1` and `HAP2` assemblies in the samplesheet, this effects naming of output files.
- Pipeline can now correctly handle a null value as a barcode value, removing the need for a dummy value in cases of data such as ONT which does not use barcodes.
- Pipeline now creates a unique list of barcodes internally, stops issues where end-user scripts collect a running list of barcodes from read files.

### `Fixed`

- Bug where in some cases the btk_run variable would not be set prior to it's use in a conditional, causing the pipeline to crash.
- Bug where btk input parameters would not be correctly set leading to incorrect runs of btk, where the wrong sample (no contamination) is used and the right sample is passed over.
- Corrected an error with output of FCSGX where it was only looking for BAM files, which are not produced by FCSGX!
- RUN_COVERAGE was using a legacy Variable name
- Bug where fcsgx was not generating output matching cobiontcheck (predecessor to ASCC), found to be caused by an incorrect threshold value.
- AUTOFILTER scripts were adapted to make the trigger values accessible to the end-user.
- Other fixes introduced a race condition for sanger-tol/btk. This is now fixed and v2 of sanger-tol/ascc will wholey remove them. #83
- Removed cpu and memory resource multipliers, it's not needed. If it crashes, it'll be for something else.

### `Bugs`

- Currently, Blobtoolkit will _not_ run if there is no autofilter output channel.
- tiara, fcsgx, autofilter must always be activated.

## v0.3.1 - Red Lamp (H1) [12/05/2025]

THIS IS STILL AN IN-DEVELOPMENT PROJECT SO THERE MAY BE BUGS.

Release 5 of sanger-tol/ascc, correcting environments and updating module structure.

### Enhancements & Fixes

- Testing for conda revealed multiple incorrect conda channels.
  - Leading to version inconsistancy between singularity/docker and conda.
- Updating module structure to be more similar to NF-Core modules.
  - All local modules now have a `environment.yml` for conda env control.
- RUN_READS_COVERAGE has stopped running due to conditionals becoming channels
  - channels cannot be compared
  - Removed channel.of(params....) to remedy this.
- Added a missing process conditional.
- KMER Analysis is now switched off in production.config
  - This primarily effects only SANGER-TOLA production
- Added new output from EXTRACT_CONTAMINANTS for parity with cobiontcheck (unreleased pre-nf_core ASCC pipeline).
  - Output will be added to `organellar_contamination_recomendations`
- Update sanger-tol/blobtoolkit to [v0.8.0 - Sprigatito](https://github.com/sanger-tol/blobtoolkit/releases/tag/0.8.0)
  - We will also be using the `miniprot` gene predictor.

### Dependencies

| Module         | Old Version | New Versions |
| -------------- | ----------- | ------------ |
| SANGER_TOL_BTK | dev         | 0.8.0        |

## v0.3.0 - Red Lamp [02/05/2025]

THIS IS STILL AN IN-DEVELOPMENT PROJECT SO THERE MAY BE BUGS.

Release 4 of sanger-tol/ascc, correcting bugs found in production testing and correcting the strucutre of the pipeline.

### Enhancements & Fixes

- Re-added the ascc.nf=.
  - This corrects an issue with version generation.
- Remove unnecessary modules which have been replaced by NF-core modules [#104](https://github.com/sanger-tol/ascc/issues/104) .
- Corrected incorrect file output
  - publishDir does not take a list of file extensions, only positive globs.
    - This has meant some modules have more patterns to match than previously.
    - Update modules which were outputting versions to the output dir unnecessarily.
- Corrected issue with channel generation being tempermental in some cases.
- Added per-process enums rather than the use of a csv list to control process execution [#107](https://github.com/sanger-tol/ascc/issues/107).
  - Significantly easier to maintain conditionals for process execution.
  - Each unique process now has a flag such as `run_{process}` which accepts a value of `genomic`,`organellar`,`both`,`off`.
  - For some processes the options are `genomic`,`both`,`off` as they arte not useful for organellar assemblies.
  - Implemented due to feedback from users.
  - This method is now significantly easier to control and understand.
  - Due to this update we are re-introducing the `genomic_only` flag.
- Update to CI and test files.
  - Includes the temporary deletion of the download_pipelines.yaml [#118](https://github.com/sanger-tol/ascc/pull/118)
  - This will require addition of databases and data to be downloaded and set up for that runner.
- Added a production profile - intended to simplify production needs in Sanger ToL [#106](https://github.com/sanger-tol/ascc/issues/106).
  - Updates configs for test and tol_assembly.
- Update modules which were requesting 100.h.
- Template has been updated to 3.2.1.
- Updated Documentation.
- Updated modules.config to output more files.
- Update blobtoolkit to be using the dev branch - essential update before btk releases. [#114](https://github.com/sanger-tol/ascc/pull/114)
- Updated FCSGX/RUNGX to version 0.5.5 [#52](https://github.com/sanger-tol/ascc/issues/52)
- Updated the conditional logic for Blobtoolkit [#123](https://github.com/sanger-tol/ascc/issues/123)
- Updates description to some fields - @prototaxites [#121](https://github.com/sanger-tol/ascc/pull/121)
- Added @prototaxites as a contributor.

### Parameters

| Old Parameter        | New Parameter             |
| -------------------- | ------------------------- |
| -                    | --genomic_only            |
| --inlcude            | REMOVED                   |
| --exclude            | REMOVED                   |
| --organellar_include | REMOVED                   |
| --organellar_exclude | REMOVED                   |
| -                    | --run_essentials          |
| -                    | --run_kmers               |
| -                    | --run_tiara               |
| -                    | --run_coverage            |
| -                    | --run_nt_blast            |
| -                    | --run_nr_diamond          |
| -                    | --run_uniprot_diamond     |
| -                    | --run_kraken              |
| -                    | --run_fcsgx               |
| -                    | --run_fcs_adaptor         |
| -                    | --run_vecscreen           |
| -                    | --run_btk_busco           |
| -                    | --run_pacbio_barcodes     |
| -                    | --run_organellar_blast    |
| -                    | --run_autofilter_assembly |
| -                    | --run_create_btk_dataset  |
| -                    | --run_merge_datasets      |

### Dependencies

| Module         | Old Version | New Versions |
| -------------- | ----------- | ------------ |
| SANGER_TOL_BTK | 0.7.1       | dev          |
| FCSGX_RUNGX    | 0.5.4       | 0.5.5        |

## v0.2.1 - Red Speaker [25/04/2025]

THIS IS STILL AN IN-DEVELOPMENT PROJECT SO THERE MAY BE BUGS.

Release 3 of sanger-tol/ascc, correcting bugs stopping use in production.

### Enhancements & Fixes

- Updating GENERATE_SAMPLESHEET module.
  - Caused by issues with test profiles.
- Updating the test profiles.
- Minor update to the base.config - 100.h was too much!
  - At sanger too many jobs were heads to the week queue, when they simply didn't need to!
- Remove unnecessary module files, they are no longer used and are replaced by NF_Core blast modules. #104
  - BLAST_V5_DATABASE -> BLAST_BLASTN
  - BLAST_MAKEBLASTDB_BARCODES -> BLAST_MAKEBLASTDB

## v0.2.0 - Red Speaker [14/04/2025]

THIS IS STILL AN IN-DEVELOPMENT PROJECT SO THERE MAY BE BUGS.

Release 2 of sanger-tol/ascc, updated with the [nf-core](https://nf-co.re/) template (v3.2).

THIS IS STILL AN IN-DEVELOPMENT PROJECT SO THERE MAY BE BUGS.

### Enhancements & Fixes

- Template Updated to 3.2 which moves us to NF-SCHEMA.
- An update to the KMER counting scripts and related processes.
- Re-organisation of .nf files into current standards.
- Updates to scripts using the ncbi_rankedlineage new_taxdump.
- Update to the GENERATE SAMPLESHEET script to convert from Python to Bash and add a counter
- Updated GET_LARGEST_SCAFFOLD to replace shell block.
- Update parse_fcsgx_result and autofilter python scripts to handle both the previous and new format of taxdump (format changed at the end of March, adding a new column for virus hierarchy).
- main.nf has been cleaned so that pipeline prep code is now found in the PIPELINE_INITIALISATION subworkflow.
- Blobtoolkit pipeline version has been updated to 0.7.1, this required changes to the input (GENERATE_SAMPLESHEET).
- Added sample_id as an input param. Currently only used in a few places.
- Adding contributors to the new format.
- Patched MAKEBLAST_DB to use it's own --out flag rather that the current mkdir and mv solution.
- Input reads now need to be specified as an array in the yaml.
- `ncbi_taxonomy_path` has been removed as the accession2taxid db is no longer necessary.
- `genomic_only` has been removed, functionality can be replicated with `--organellar_exclude ALL`.

### Parameters

| Old Parameter        | New Parameter |
| -------------------- | ------------- |
| -                    | --sample_id   |
| --ncbi_taxonomy_path | REMOVED       |
| --genomic_only       | REMOVED       |

### Dependencies

| Module                        | Old Version                                                            | New Versions                                                           |
| ----------------------------- | ---------------------------------------------------------------------- | ---------------------------------------------------------------------- |
| PARSE_FCSGX_RESULT            | python:3.9, parse_fcsgx_result.py:1.0.0                                | python:3.9, parse_fcsgx_result.py:1.0.1                                |
| AUTOFILTER_AND_CHECK_ASSEMBLY | python:3.9, autofilter.py:1.0.0, abnormal_contamination_check.py:1.0.0 | python:3.9, autofilter.py:1.0.0, abnormal_contamination_check.py:1.0.1 |
| BLASTN                        | 2.14.0                                                                 | 2.16.0                                                                 |
| MAKEBLASTDB                   | 2.14.0                                                                 | 2.16.0                                                                 |
| GENERATE_SAMPLESHEET          | 1.0.0                                                                  | 1.1.0                                                                  |

## v0.1.0 - Red Book [14/02/2025]

Initial release of sanger-tol/ascc, created with the [nf-core](https://nf-co.re/) template.

THIS IS STILL AN IN DEVELOPMENT PROJECT SO THERE MAY BE BUGS.

The intention of this pipeline is to succeed the currently in production CobiontCheck (A version of ASCC which is not NF-core Compliant and not publically available).

### Enhancements & Fixes

- Subworkflow to generate NT_BLAST results
- Subworkflow to generate TIARA hits results.
- Subworkflow to generate a .genome file.
- Subworkflow to generate KMER profiles.
- Subworkflow to generate BLAST results of Organellar seqeuence against genomic.
- Subworkflow to generate MINIMAP2 mapping of reads.
- Subworkflow to generate COVERAGE data across the genome.
- Subworkflow to generate DIAMOND blast results.
- Subworkflow to generate PACBIO_BARCODE_CHECK results using BLAST.
- Subworkflow to generate FCS_ADAPTOR results.
- Subworkflow to generate FCS_GX results.
- Subworkflow to generate VECSCREEN results.
- Subworkflow to remove trailing sequences of N's.
- Subpipeline to run the SANGER_TOL_BTK pipeline.
- Module to automatically filter the genome.
- Modules for generating and merging BTK dataset folders.
- Module for generating samplesheets required for SANGER_TOL_BTK.
- Added preliminary paper for JOSS.
- Added documentation.
- Updated config files.
- Updated the CICD config files.

### Parameters

| Old Parameter | New Parameter                      |
| ------------- | ---------------------------------- |
| -             | --include                          |
| -             | --exclude                          |
| -             | --organellar_include               |
| -             | --organellar_include               |
| -             | --genomic_only                     |
| -             | --reads_path                       |
| -             | --reads_type                       |
| -             | --pacbio_barcode_file              |
| -             | --pacbio_barcode_names             |
| -             | --scientific_name                  |
| -             | --taxid                            |
| -             | --kmer_length                      |
| -             | --dimensionality_reduction_methods |
| -             | --nt_database_path                 |
| -             | --nt_database_prefix               |
| -             | --nt_kraken_database_path          |
| -             | --ncbi_accession_ids_folder        |
| -             | --ncbi_taxonomy_path               |
| -             | --ncbi_ranked_lineage_path         |
| -             | --busco_lineages_folder            |
| -             | --busco_lineages                   |
| -             | --fcs_gx_database_path             |
| -             | --vecscreen_database_path          |
| -             | --diamond_uniprot_database_path    |
| -             | --diamond_nr_database_path         |
| -             | --seqkit_sliding                   |
| -             | --seqkit_window                    |
| -             | --btk_yaml                         |
| -             | --n_neighbours                     |
| -             | --btk_busco_run_mode               |

### Dependencies

| Module                                  | Old Version | New Versions                                                                                                                      |
| --------------------------------------- | ----------- | --------------------------------------------------------------------------------------------------------------------------------- |
| BLASTN                                  | -           | 2.15.0--pl5321h6f7f691_1                                                                                                          |
| MAKEBLASTDB                             | -           | 2.15.0--pl5321h6f7f691_1                                                                                                          |
| CUSTOM_GETCHROMSIZES                    | -           | samtools:1.21--h50ea8bc_0                                                                                                         |
| DIAMOND_BLASTX                          | -           | 2.1.8--h43eeafb_0                                                                                                                 |
| FCS_FCSADAPTOR                          | -           | 0.5.0                                                                                                                             |
| FCSGX_RUNGX                             | -           | 0.5.4--h4ac6f70_1                                                                                                                 |
| GNU_SORT                                | -           | 9.3                                                                                                                               |
| GUNZIP                                  | -           | ubuntu:22.04                                                                                                                      |
| KRAKEN2_KRAKEN2                         | -           | kraken2:2.1.3,pigz:2.8                                                                                                            |
| MINIMAP2_ALIGN                          | -           | minimap2:2.28--he4a0461_0,samtools=1.20                                                                                           |
| MINIMAP2_INDEX                          | -           | 2.28--he4a0461_0                                                                                                                  |
| NCBITOOLS_VECSCREEN                     | -           | ncbi-tools-bin:6.1.20170106-6-deb_cv2                                                                                             |
| SAMTOOLS\_\*                            | -           | 1.21--h50ea8bc_0                                                                                                                  |
| SEQKIT_SLIDING                          | -           | 2.8.1--h9ee0642_0                                                                                                                 |
| TIARA                                   | -           | 1.0.3                                                                                                                             |
| MAKEBLASTDB_PACBIO_BARCODES             | -           | 2.15.0--pl5321h6f7f691_1                                                                                                          |
| SANGER_TOL_BTK                          | -           | 0.6.0 Bellsprout                                                                                                                  |
| ASCC_MERGE_TABLES                       | -           | pandas:1.5.2, python:3.10, ascc_merge_tables.py:1.0.0                                                                             |
| AUTOFILTER_AND_CHECK_ASSEMBLY           | -           | python:3.9, autofilter.py:1.0.0, abnormal_contamination_check.py:1.0.0                                                            |
| BLAST_CHUNK_TO_FULL                     | -           | python:3.9, blast_hit_chunk_coords_to_full_coords.py:1.0.0                                                                        |
| BLAST_GET_TOP_HITS                      | -           | pandas:1.5.2, python:3.10, blast_get_top_hits.py:1.0.0                                                                            |
| CHECK_BARCODE                           | -           | python:3.9. pacbio_barcode_check.py:1.0.0                                                                                         |
| CHUNK_ASSEMBLY_FOR_VECSCREEN            | -           | biopython:1.81, chunk_assembly_for_vecscreen.py:1.0.0                                                                             |
| CONVERT_TO_HITS_FILE                    | -           | python:3.9, convert_to_hits.py:1.0.0                                                                                              |
| CREATE_BTK_DATASET                      | -           | blobtoolkit:4.3.9, python:3.9, create_btk_dataset.py:2.0.0                                                                        |
| EXTRACT_CONTAMINANTS                    | -           | python:3.9, biopython:1.78, pybedtools:0.9.0, extract_contaminants_by_type.py:1.0.0                                               |
| FILTER_BARCODE                          | -           | biopython:1.78, python:3.9, filter_barcode_blast_results.py:1.0.0                                                                 |
| FILTER_COMMENTS                         | -           | coreutils:9.1                                                                                                                     |
| FILTER_FASTA                            | -           | python:3.9, sanitise_input_fasta_file.py:1.2.0, filter_fasta_by_length.py:1.0.0                                                   |
| FILTER_VECSCREEN_RESULTS                | -           | python:3.9, VSlistTo1HitPerLine.py:1.0.0                                                                                          |
| REFORMAT_DIAMOND_OUTFMT6                | -           | python:3.9, reformat_diamond_outfmt6.py:1.0.0                                                                                     |
| GC_CONTENT                              | -           | python:3.9, gc_content.py:1.0.0                                                                                                   |
| GENERATE_SAMPLESHEET                    | -           | python:3.9, generate_samplesheet.py:1.0.0                                                                                         |
| GET_KMER_COUNTS                         | -           | python:3.9, kcounter:0.1.1, get_kmers_counts.py:1.0.0, general_purpose_functions.py:1.0.0                                         |
| GET_LARGEST_SCAFF                       | -           | coreutils:9.1                                                                                                                     |
| GET_LINEAGE_FOR_KRAKEN                  | -           | pandas:1.5.2, python:3.9, general_purpose_functions.py:1.0.0, get_lineage_for_kraken_results.py:1.0.0                             |
| GET_LINEAGE_FOR_TOP                     | -           | python:3.9, get_lineage_for_top.py:1.0.0                                                                                          |
| KMER_COUNT_DIM_REDUCTION_COMBINE_CSV    | -           | pandas:1.5.2, python:3.9, kmer_count_dim_reduction_combine_csv.py:1.0.0                                                           |
| KMER_COUNT_DIM_REDUCTION                | -           | python:3.9, pandas:2.2.1, tensorlflow:2.15.0, scikit-learn:1.4.1, umap:0.5.5, matplotlib:3.8.0, kmer_count_dim_reduction.py:1.0.0 |
| MERGE_BTK_DATASETS                      | -           | blobtoolkit:4.3.9, merge_btk_datasets.py:2.0.0                                                                                    |
| ORGANELLE_CONTAMINATION_RECOMMENDATIONS | -           | python:3.9, organelle_contamination_recommendation.py:1.0.0                                                                       |
| PARSE_FCSGX_RESULT                      | -           | python:3.9, parse_fcsgx_result.py:1.0.0                                                                                           |
| REFORMAT_FULL_OUTFMT6                   | -           | python:3.9, reformat_blast_outfmt6.py:1.0.0                                                                                       |
| SAMTOOLS_DEPTH_AVERAGE_COVERAGE         | -           | python:3.9, samtools_depth_average_coverage.py:1.0.0                                                                              |
| SED_SED                                 | -           | ubuntu:20.04                                                                                                                      |
| SUMMARISE_VECSCREEN_OUTPUT              | -           | python:3.9, summarise_vecscreen_output.py:1.0.0                                                                                   |
| TRAILINGNS                              | -           | biopython:1.81, python:3.9, trim_Ns.py:1.0.0                                                                                      |
| VALIDATE_TAXID                          | -           | python:3.9, find_taxid_in_taxdump.py:1.0.0                                                                                        |
