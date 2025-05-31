/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ASCC GENOMIC REPORTING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CREATE_BTK_DATASET                            } from '../modules/local/blobtoolkit/create_dataset/main'
include { MERGE_BTK_DATASETS                            } from '../modules/local/blobtoolkit/merge_dataset/main'
include { ASCC_MERGE_TABLES                             } from '../modules/local/ascc/merge_tables/main'
include { AUTOFILTER_AND_CHECK_ASSEMBLY                 } from '../modules/local/autofilter/autofilter/main'
include { SANGER_TOL_BTK                                } from '../modules/local/sanger-tol/btk/main'
include { GENERATE_SAMPLESHEET                          } from '../modules/local/blobtoolkit/generate_samplesheet/main'
include { NEXTFLOW_RUN as SANGER_TOL_BTK_CASCADE        } from '../modules/local/run/main'

workflow ASCC_GENOMIC_REPORTING {
    take:
    reference_tuple_from_GG      // Channel: reference tuple from GENERATE_GENOME
    ej_dot_genome                // Channel: dot_genome from ESSENTIAL_JOBS
    ej_gc_coverage               // Channel: gc_content_txt from ESSENTIAL_JOBS
    include_workflow_steps       // List of included workflow steps
    exclude_workflow_steps       // List of excluded workflow steps
    kmers                        // Channel: from analysis workflow
    tiara                        // Channel: from analysis workflow
    nt_blast                     // Channel: from analysis workflow
    blast_lineage                // Channel: from analysis workflow
    btk_format                   // Channel: from analysis workflow
    nr_full                      // Channel: from analysis workflow
    nr_hits                      // Channel: from analysis workflow
    un_full                      // Channel: from analysis workflow
    un_hits                      // Channel: from analysis workflow
    fcsgx                        // Channel: from analysis workflow
    coverage                     // Channel: from analysis workflow
    bam                          // Channel: from analysis workflow
    vecscreen                    // Channel: from analysis workflow
    kraken1                      // Channel: from analysis workflow
    kraken2                      // Channel: from analysis workflow
    kraken3                      // Channel: from analysis workflow
    fcsadapt                     // Channel: from analysis workflow
    scientific_name              // val(name)
    btk_busco_run_mode          // val(mode)
    trailing_ns_report          // Channel: from essential jobs
    filter_fasta_sanitation_log // Channel: from essential jobs
    filter_fasta_length_filtering_log // Channel: from essential jobs
    pacbio_barcode_check_filtered // Channel: from analysis workflow
    fcsgx_report                // Channel: from analysis workflow
    fcsgx_taxonomy_report       // Channel: from analysis workflow
    kmers_results               // Channel: from analysis workflow

    main:
    ch_versions = Channel.empty()

    // Reporting workflows will be extracted from the monolithic file
    // This is a placeholder structure - content will be added in next step

    emit:
    versions = ch_versions
}
