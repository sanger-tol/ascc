#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    sanger-tol/ascc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/sanger-tol/ascc
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ASCC_GENOMIC                                      } from './workflows/ascc_genomic'
include { ASCC_ORGANELLAR                                   } from './workflows/ascc_organellar'

include { PIPELINE_INITIALISATION                           } from './subworkflows/local/utils_nfcore_ascc_pipeline'
include { PIPELINE_COMPLETION                               } from './subworkflows/local/utils_nfcore_ascc_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow SANGERTOL_ASCC_GENOMIC {

    take:
    samplesheet // channel: samplesheet read in from --input
    organelles
    include_steps
    exclude_steps
    fcs
    read_files
    scientific_name
    pacbio_db

    main:

    //
    // WORKFLOW: Run pipeline
    //
    ASCC_GENOMIC (
        samplesheet,
        organelles,
        include_steps,
        exclude_steps,
        fcs,
        read_files,
        scientific_name,
        pacbio_db
    )
}

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow SANGERTOL_ASCC_ORGANELLAR {

    take:
    samplesheet // channel: samplesheet read in from --input
    include_steps
    exclude_steps
    fcs
    reads
    scientific_name
    pacbio_db

    main:

    //
    // WORKFLOW: Run pipeline
    //
    ASCC_ORGANELLAR (
        samplesheet,
        include_steps,
        exclude_steps,
        fcs,
        reads,
        scientific_name,
        pacbio_db
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )


    //
    // WORKFLOW: Run main workflow for GENOMIC samples
    //
    // SANGERTOL_ASCC_GENOMIC (
    //     PIPELINE_INITIALISATION.out.main_genomes,
    //     PIPELINE_INITIALISATION.out.organellar_genomes,
    //     PIPELINE_INITIALISATION.out.include_steps,
    //     PIPELINE_INITIALISATION.out.exclude_steps,
    //     PIPELINE_INITIALISATION.out.fcs_gx_database,
    //     PIPELINE_INITIALISATION.out.collected_reads,
    //     Channel.of(params.scientific_name),
    //     PIPELINE_INITIALISATION.out.pacbio_db,
    // )


    //
    // WORKFLOW: Run main workflow for ORGANELLAR samples
    //
    include_workflow_steps  = params.organellar_include ? params.organellar_include.split(",") : "ALL"
    exclude_workflow_steps  = params.organellar_exclude ? params.organellar_exclude.split(",") : "NONE"

    if ( exclude_workflow_steps.contains('ALL') || include_workflow_steps.contains('NONE') ) {
        log.warn "ORGANELLAR SUBWORKFLOW: HAS BEEN SKIPPED: $exclude_workflow_steps"

    } else {
        SANGERTOL_ASCC_ORGANELLAR (
            PIPELINE_INITIALISATION.out.organellar_genomes,
            PIPELINE_INITIALISATION.out.organellar_include,
            PIPELINE_INITIALISATION.out.organellar_exclude,
            PIPELINE_INITIALISATION.out.fcs_gx_database,
            PIPELINE_INITIALISATION.out.collected_reads,
            Channel.of(params.scientific_name),
            PIPELINE_INITIALISATION.out.pacbio_db,
        )
    }


    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        []
    )
}

process MAIN_WORKFLOW_GrabFiles {
    tag "Grab PacBio Data"
    executor 'local'

    input:
    path("in")

    output:
    path("in/*.{fa,fasta,fna}.{gz}")

    "true"
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
