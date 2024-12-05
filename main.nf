#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    sanger-tol/ascc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/sanger-tol/ascc
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VALIDATE_TAXID            } from './modules/local/validate_taxid'
include { COPY_FCS_GX               } from './modules/local/copy_fcs_tmp'

include { ASCC_GENOMIC              } from './workflows/ascc_genomic'
include { ASCC_ORGANELLAR           } from './workflows/ascc_organellar'

include { PIPELINE_INITIALISATION   } from './subworkflows/local/utils_nfcore_ascc_pipeline'
include { PIPELINE_COMPLETION       } from './subworkflows/local/utils_nfcore_ascc_pipeline'

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
    validate_versions
    include_steps
    exclude_steps
    fcs

    main:

    //
    // WORKFLOW: Run pipeline
    //
    ASCC_GENOMIC (
        samplesheet,
        organelles,
        validate_versions,
        include_steps,
        exclude_steps,
        fcs
    )
}

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow SANGERTOL_ASCC_ORGANELLAR {

    take:
    samplesheet // channel: samplesheet read in from --input
    validate_versions
    include_steps
    exclude_steps
    fcs

    main:

    //
    // WORKFLOW: Run pipeline
    //
    ASCC_ORGANELLAR (
        samplesheet,
        validate_versions,
        include_steps,
        exclude_steps,
        fcs
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    //
    // WORKFLOW: THIS WHOLE THING SHOULD BE IN A SUBWORKFLOW REALY
    //

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // LOGIC: FILTER THE INPUT BASED ON THE assembly_type VALUE IN THE META
    //          DEPENDING ON THIS VALUE THE PIPELINE WILL NEED TO BE DIFFERENT
    //
    PIPELINE_INITIALISATION.out.samplesheet
        .branch{
            organellar_genome: it[0].assembly_type == "MITO" || it[0].assembly_type == "PLASTID"
            sample_genome: it[0].assembly_type  == "PRIMARY" || it[0].assembly_type  == "HAPLO"
            error: true
        }
        .set { branched_assemblies }


    //
    // MODULE: ENSURE THAT THE TAXID FOR THE INPUT GENOME IS INDEED IN THE TAXDUMP
    //
    VALIDATE_TAXID(
        params.taxid,
        params.ncbi_taxonomy_path
    )


    //
    // TEMP MODULE: COPY FCS_GX_DB INTO TMP
    //
    include_workflow_steps  = params.include ? params.include.split(",") : "ALL"
    VALID = ["ALL", "fcs-gx"]

    if ( params.move_fcs && VALID.any{ include_workflow_steps.contains( it ) } ) {
        COPY_FCS_GX (
                params.fcs_gx_database_path,
        )

        fcs_gx_database_path = COPY_FCS_GX.out.fcsdb_path
        delete_fcs = true
    } else {
        fcs_gx_database_path = Channel.of(params.fcs_gx_database_path)
        delete_fcs = false
    }


    //
    // WORKFLOW: Run main workflow for GENOMIC samples
    //
    SANGERTOL_ASCC_GENOMIC (
        branched_assemblies.sample_genome,
        branched_assemblies.organellar_genome,
        VALIDATE_TAXID.out.versions,
        params.include,
        params.exclude,
        fcs_gx_database_path,
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
    SANGERTOL_ASCC_ORGANELLAR (
        branched_assemblies.organellar_genome,
        VALIDATE_TAXID.out.versions,
        organellar_include,
        organellar_exclude,
        fcs_gx_database_path
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
