/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ESSENTIAL_JOBS                                } from '../subworkflows/local/essential_jobs/main'
include { ASCC_GENOMIC_ANALYSIS                         } from './ascc_genomic_analysis'
include { ASCC_GENOMIC_REPORTING                        } from './ascc_genomic_reporting'

include { paramsSummaryMultiqc                          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                        } from '../subworkflows/local/utils_nfcore_ascc_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASCC_GENOMIC {

    take:
    ch_samplesheet          // channel: samplesheet read in from --input
    organellar_genomes      // channel: tuple(meta, reference)
    fcs_db                  // [path(path)]
    reads
    scientific_name         // val(name)
    pacbio_database         // tuple [[meta.id], pacbio_database]
    ncbi_taxonomy_path
    ncbi_ranked_lineage_path
    nt_database_path
    diamond_nr_db_path
    diamond_uniprot_db_path
    taxid                   // NEW PARAMETER from dev branch
    nt_kraken_db_path       // NEW PARAMETER from dev branch
    vecscreen_database_path
    reads_path
    reads_layout
    reads_type
    btk_lineages
    btk_lineages_path

    main:
    ch_versions = Channel.empty()

    //
    // LOGIC: USING DEV BRANCH run_* PARAMETER SYSTEM
    //
    // log.info "GENOMIC RUN -- Using run_* parameter system from dev branch"

    //
    // LOGIC: CREATE btk_busco_run_mode VALUE
    //
    btk_busco_run_mode = params.btk_busco_run_mode ?: "conditional"

    //
    // LOGIC: PRETTY NOTIFICATION OF FILES AT STAGE
    //
    ch_samplesheet
        .map { meta, sample ->
            log.info "GENOMIC WORKFLOW:\n\t-- $meta\n\t-- $sample"
        }

    //
    // SUBWORKFLOW: RUNS FILTER_FASTA, GENERATE .GENOME, CALCS GC_CONTENT AND FINDS RUNS OF N's
    //                  THIS SHOULD NOT RUN ONLY WHEN SPECIFICALLY REQUESTED
    //
    if (params.run_essentials == "both" ||
        (params.run_essentials == "genomic" && params.genomic_only) ||
        (params.run_essentials == "organellar" && !params.genomic_only)) {

        ESSENTIAL_JOBS(
            ch_samplesheet
        )
        ch_versions = ch_versions.mix(ESSENTIAL_JOBS.out.versions)

        reference_tuple_from_GG = ESSENTIAL_JOBS.out.reference_tuple_from_GG
        ej_dot_genome           = ESSENTIAL_JOBS.out.dot_genome
        ej_gc_coverage          = ESSENTIAL_JOBS.out.gc_content_txt
        reference_tuple_w_seqkt = ESSENTIAL_JOBS.out.reference_with_seqkit
        trailing_ns_report      = ESSENTIAL_JOBS.out.trailing_ns_report
        filter_fasta_sanitation_log = ESSENTIAL_JOBS.out.filter_fasta_sanitation_log
        filter_fasta_length_filtering_log = ESSENTIAL_JOBS.out.filter_fasta_length_filtering_log

    } else {
        log.warn("MAKE SURE YOU ARE AWARE YOU ARE SKIPPING ESSENTIAL JOBS, THIS INCLUDES BREAKING SCAFFOLDS OVER 1.9GB, FILTERING N\'s AND GC CONTENT REPORT (THIS WILL BREAK OTHER PROCESSES AND SHOULD ONLY BE RUN WITH `--run_essentials genomic`)")

        reference_tuple_from_GG = ch_samplesheet // This is the reference genome input channel
        ej_dot_genome = Channel.of([[],[]])
        ej_gc_coverage = Channel.of([[],[]])
        reference_tuple_w_seqkt = Channel.of([[],[]])
        trailing_ns_report = Channel.of([[],[]])
        filter_fasta_sanitation_log = Channel.of([[],[]])
        filter_fasta_length_filtering_log = Channel.of([[],[]])
    }

    //
    // SUBWORKFLOW: RUN ANALYSIS STEPS
    //
    ASCC_GENOMIC_ANALYSIS(
        reference_tuple_from_GG,
        organellar_genomes,
        fcs_db,
        reads,
        pacbio_database,
        ncbi_taxonomy_path,
        ncbi_ranked_lineage_path,
        nt_database_path,
        diamond_nr_db_path,
        diamond_uniprot_db_path,
        taxid,
        nt_kraken_db_path
    )
    ch_versions = ch_versions.mix(ASCC_GENOMIC_ANALYSIS.out.versions)

    // Initialize channels for PACBIO_BARCODE_CHECK and RUN_FCSGX outputs
    pacbio_barcode_check_filtered = Channel.of([[],[]])
    fcsgx_report = Channel.of([[],[]])
    fcsgx_taxonomy_report = Channel.of([[],[]])

    // Get PACBIO_BARCODE_CHECK outputs if the workflow was run
    if (params.run_pacbio_barcodes == "both" ||
        (params.run_pacbio_barcodes == "genomic" && params.genomic_only) ||
        (params.run_pacbio_barcodes == "organellar" && !params.genomic_only)) {
        pacbio_barcode_check_filtered = ASCC_GENOMIC_ANALYSIS.out.barcode_check_filtered
    }

    // Get RUN_FCSGX outputs if the workflow was run
    if (params.run_fcsgx == "both" ||
        (params.run_fcsgx == "genomic" && params.genomic_only) ||
        (params.run_fcsgx == "organellar" && !params.genomic_only)) {
        fcsgx_report = ASCC_GENOMIC_ANALYSIS.out.fcsgx_report_txt
        fcsgx_taxonomy_report = ASCC_GENOMIC_ANALYSIS.out.fcsgx_taxonomy_rpt
    }

    //
    // SUBWORKFLOW: RUN REPORTING STEPS
    //
    ASCC_GENOMIC_REPORTING(
        reference_tuple_from_GG,
        ej_dot_genome,
        ej_gc_coverage,
        ASCC_GENOMIC_ANALYSIS.out.kmers,
        ASCC_GENOMIC_ANALYSIS.out.tiara,
        ASCC_GENOMIC_ANALYSIS.out.nt_blast,
        ASCC_GENOMIC_ANALYSIS.out.blast_lineage,
        ASCC_GENOMIC_ANALYSIS.out.btk_format,
        ASCC_GENOMIC_ANALYSIS.out.nr_full,
        ASCC_GENOMIC_ANALYSIS.out.nr_hits,
        ASCC_GENOMIC_ANALYSIS.out.un_full,
        ASCC_GENOMIC_ANALYSIS.out.un_hits,
        ASCC_GENOMIC_ANALYSIS.out.fcsgx,
        ASCC_GENOMIC_ANALYSIS.out.coverage,
        ASCC_GENOMIC_ANALYSIS.out.bam,
        ASCC_GENOMIC_ANALYSIS.out.vecscreen,
        ASCC_GENOMIC_ANALYSIS.out.kraken1,
        ASCC_GENOMIC_ANALYSIS.out.kraken2,
        ASCC_GENOMIC_ANALYSIS.out.kraken3,
        ASCC_GENOMIC_ANALYSIS.out.fcsadapt,
        scientific_name,
        btk_busco_run_mode,
        trailing_ns_report,
        filter_fasta_sanitation_log,
        filter_fasta_length_filtering_log,
        pacbio_barcode_check_filtered,
        fcsgx_report,
        fcsgx_taxonomy_report,
        ASCC_GENOMIC_ANALYSIS.out.kmers_results
    )
    ch_versions = ch_versions.mix(ASCC_GENOMIC_REPORTING.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    emit:
    versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
