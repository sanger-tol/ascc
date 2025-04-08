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
    reads
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
        reads,
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


    include_workflow_steps  = params.include ? params.include.split(",") : "ALL"
    exclude_workflow_steps  = params.exclude ? params.exclude.split(",") : "NONE"


    //
    // LOGIC: GETS PACBIO READ PATHS FROM READS_PATH IF (COVERAGE OR BTK SUBWORKFLOW IS ACTIVE) OR ALL
    //
    if (
        (
            (include_workflow_steps.contains('coverage') && !exclude_workflow_steps.contains("coverage")) ||
            (include_workflow_steps.contains('btk_busco') && !exclude_workflow_steps.contains("btk_busco"))
        ) || (
            include_workflow_steps.contains('ALL') && !exclude_workflow_steps.contains("btk_busco") && !exclude_workflow_steps.contains("coverage")
        ) || (
            include_workflow_steps.contains('ALL')
        )
    ) {
        ch_grabbed_reads_path       = Channel.of(params.reads_path).collect()
    } else {
        ch_grabbed_reads_path       = []
    }


    //
    // WORKFLOW: Run main workflow for GENOMIC samples
    //
    // TODO: THIS WOULD HAVE BEEN SIMPLER TO FIX BY COMBINING THE ORGANELLAR GENOMES TO GENOMIC!!!
    SANGERTOL_ASCC_GENOMIC (
        PIPELINE_INITIALISATION.out.main_genomes,
        PIPELINE_INITIALISATION.out.organellar_genomes,
        params.include,
        params.exclude,
        PIPELINE_INITIALISATION.out.fcs_gx_database,
        ch_grabbed_reads_path,
        Channel.of(params.scientific_name),
        PIPELINE_INITIALISATION.out.pacbio_db,
    )


    //
    // LOGIC: IT NOW MAKES SENSE TO BE USING SPECIFIC FLAGS FOR THE ORGANELLAR WORKFLOW
    //          IF THEY ARE NOT SPECIFIED THEN THE USE OF THE GENOMIC ONES WILL SUFFICE.
    //
    if ( !params.organellar_include && params.include ) {
        println "Using GENOMIC specific include/exclude flags (make sure you are supposed to be!)"
        organellar_include = params.include
    } else {
        println "Using ORGANELLE specific include/exclude flags"
        organellar_include = params.organellar_include
    }

    if ( !params.organellar_exclude && params.exclude ) {
        println "Using GENOMIC specific include/exclude flags (make sure you are supposed to be!)"
        organellar_exclude = params.exclude
    } else {
        println "Using ORGANELLE specific include/exclude flags"
        organellar_exclude = params.organellar_exclude
    }


    //
    // WORKFLOW: Run main workflow for ORGANELLAR samples
    //

    if ( !params.genomic_only ) {

        SANGERTOL_ASCC_ORGANELLAR (
            branched_assemblies.organellar_genome,
            organellar_include,
            organellar_exclude,
            PIPELINE_INITIALISATION.out.fcs_gx_database,
            ch_grabbed_reads_path,
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
