# sanger-tol/ascc: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.3.1 - Red Lamp (H1) [12/05/2025]

Release 5 of sanger-tol/ascc, correcting environments and updating module structure.

### Enhancements & Fixes

- Testing for conda revealed a few incorrect conda channels.
- Updating module structure to be more similar to NF-Core modules.
- RUN_READS_COVERAGE has stopped running due to conditionals becoming channels
  - channels cannot be compared
  - Removed channel.of(params....) to remedy this.

## v0.3.0 - Red Lamp [02/05/2025]

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
