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

    main:
    ch_versions = Channel.empty()

    //
    // LOGIC: CONTROL OF THE INCLUDE AND EXCLUDE FLAGS
    //      TODO: THESE SHOULD CREATE A SET OF INCLUDE - EXCLUDE
    //      TODO: YES THIS IS DUPLICATED FROM PIPELINE INIT,
    //              HOWEVER THAT CONVERTED THE VALUES INTO A CHANNEL WHICH ISN'T THE EASIEST THING TO THEN PARSE OUT
    include_workflow_steps  = params.include ? params.include.split(",") : "ALL"
    exclude_workflow_steps  = params.exclude ? params.exclude.split(",") : "NONE"

    full_list               = [
        "essentials", "kmers", "tiara", "coverage", "nt_blast", "nr_diamond",
        "uniprot_diamond", "kraken", "fcs-gx", "fcs-adaptor", "vecscreen", "btk_busco",
        "pacbio_barcodes", "organellar_blast", "autofilter_assembly", "create_btk_dataset",
        "merge", "html_report", "ALL", "NONE"
    ]

    if (!full_list.containsAll(include_workflow_steps) && !full_list.containsAll(exclude_workflow_steps)) {
        exit 1, "There is an extra argument given on Command Line: \n Check contents of: $include_workflow_steps\nAnd $exclude_workflow_steps\nMaster list is: $full_list"
    }

    log.info "GENOMIC RUN -- INCLUDE STEPS INC.: $include_workflow_steps"
    log.info "GENOMIC RUN -- EXCLUDE STEPS INC.: $exclude_workflow_steps"

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
    if ( !exclude_workflow_steps.contains("essentials") && (include_workflow_steps.contains("ALL") || include_workflow_steps.contains("essentials")) ) {

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
        log.warn("MAKE SURE YOU ARE AWARE YOU ARE SKIPPING ESSENTIAL JOBS, THIS INCLUDES BREAKING SCAFFOLDS OVER 1.9GB, FILTERING N\'s AND GC CONTENT REPORT (THIS WILL BREAK OTHER PROCESSES AND SHOULD ONLY BE RUN WITH `--include essentials`)")

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
        include_workflow_steps,
        exclude_workflow_steps,
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
    if ((include_workflow_steps.contains('pacbio_barcodes') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("pacbio_barcodes")) {
        pacbio_barcode_check_filtered = ASCC_GENOMIC_ANALYSIS.out.barcode_check_filtered
    }

    // Get RUN_FCSGX outputs if the workflow was run
    if ((include_workflow_steps.contains('fcs-gx') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("fcs-gx")) {
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
        include_workflow_steps,
        exclude_workflow_steps,
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
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
