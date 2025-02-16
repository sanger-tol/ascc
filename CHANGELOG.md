# sanger-tol/ascc: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
