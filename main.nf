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

include { VALIDATE_TAXID as MAIN_WORKFLOW_VALIDATE_TAXID    } from './modules/local/validate_taxid'
include { GUNZIP as MAIN_WORKFLOW_GUNZIP                    } from './modules/nf-core/gunzip/main'
include { CHECK_NT_BLAST_TAXONOMY                           } from './modules/local/check_nt_blast_taxonomy/main'

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
    validate_versions
    include_steps
    exclude_steps
    fcs
    reads

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
        fcs,
        reads
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
    reads

    main:

    //
    // WORKFLOW: Run pipeline
    //
    ASCC_ORGANELLAR (
        samplesheet,
        validate_versions,
        include_steps,
        exclude_steps,
        fcs,
        reads
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
    ch_versions = Channel.empty()


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

    fcs_gx_database_path = Channel.of(params.fcs_gx_database_path)

    //
    // LOGIC: GUNZIP INPUT DATA IF GZIPPED, OTHERWISE PASS
    //

    PIPELINE_INITIALISATION.out.samplesheet
        .branch { meta, file ->
            zipped: file.name.endsWith('.gz')
            unzipped: !file.name.endsWith('.gz')
        }
        .set {ch_input}


    //
    // MODULE: UNZIP INPUTS IF NEEDED
    //
    MAIN_WORKFLOW_GUNZIP (
        ch_input.zipped
    )

    // TODO: NOT THE RIGHT PLACE FOR THIS
    // SHOULD THIS BE MOVED INTO A SUBWORKFLOW AND EXECUTED INSIDE THE SUBWORKFLOWS
    ch_versions = ch_versions.mix(MAIN_WORKFLOW_GUNZIP.out.versions)


    //
    // LOGIC: MIX CHANELS WHICH MAY OR MAY NOT BE EMPTY INTO A SINGLE QUEUE CHANNEL
    //
    unzipped_input = Channel.empty()

    unzipped_input
        .mix(ch_input.unzipped, MAIN_WORKFLOW_GUNZIP.out.gunzip)
        .set { standardised_unzipped_input }


    //
    // LOGIC: FILTER THE INPUT BASED ON THE assembly_type VALUE IN THE META
    //          DEPENDING ON THIS VALUE THE PIPELINE WILL NEED TO BE DIFFERENT
    //
    standardised_unzipped_input
        .branch{
            organellar_genome: it[0].assembly_type == "MITO" || it[0].assembly_type == "PLASTID"
            sample_genome: it[0].assembly_type  == "PRIMARY" || it[0].assembly_type  == "HAPLO"
            error: true
        }
        .set { branched_assemblies }


    include_workflow_steps  = params.include ? params.include.split(",") : "ALL"
    exclude_workflow_steps  = params.exclude ? params.exclude.split(",") : "NONE"


    //
    // MODULE: ENSURE THAT THE TAXID FOR THE INPUT GENOME IS INDEED IN THE TAXDUMP
    //
    MAIN_WORKFLOW_VALIDATE_TAXID(
        params.taxid,
        params.ncbi_taxonomy_path
    )

    //
    // LOGIC: PARSE WORKFLOW STEPS FOR RUN
    //
    
    include_workflow_steps_genomic = params.include ? params.include.split(",") : ["ALL"]
    exclude_workflow_steps_genomic = params.exclude ? params.exclude.split(",") : ["NONE"]

    include_workflow_steps_organellar = params.organellar_include ? params.organellar_include.split(",") : include_workflow_steps_genomic
    exclude_workflow_steps_organellar = params.organellar_exclude ? params.organellar_exclude.split(",") : exclude_workflow_steps_genomic

    //
    // LOGIC: CHECK IF NT BLAST IS INCLUDED IN EITHER GENOMIC OR ORGANELLAR WORKFLOW
    //

    run_nt_blast_genomic = (include_workflow_steps_genomic.contains('nt_blast') || include_workflow_steps_genomic.contains('ALL')) && !exclude_workflow_steps_genomic.contains("nt_blast")
    run_nt_blast_organellar = (include_workflow_steps_organellar.contains('nt_blast') || include_workflow_steps_organellar.contains('ALL')) && !exclude_workflow_steps_organellar.contains("nt_blast")

    //
    // MODULE: CHECK IF NT BLAST DATABASE HAS TAXONOMY INCLUDED (ONLY IF NT BLAST IS INCLUDED)
    // This check is specifically for the nt BLAST database used in the EXTRACT_NT_BLAST subworkflow,
    // not for other BLAST databases used elsewhere in the pipeline (VecScreen, PacBio barcodes check, etc.)
    //
    if (run_nt_blast_genomic || run_nt_blast_organellar) {
        CHECK_NT_BLAST_TAXONOMY(
            params.nt_database_path
        )
        ch_versions = ch_versions.mix(CHECK_NT_BLAST_TAXONOMY.out.versions)
        
        // Check the result and fail if needed
        CHECK_NT_BLAST_TAXONOMY.out.status
            .map { it.trim() }  // Trim any whitespace
            .subscribe { status ->
                if (status == "nt_database_taxonomy_files_not_found") {
                    log.error "NT BLAST database taxonomy check failed"
                    exit 1, "The NT BLAST database does not have taxonomy included. Please see the error message above for details."
                }
            }
    }

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
            include_workflow_steps.contains('ALL') && params.profile_name == 'test'
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
        branched_assemblies.sample_genome,
        branched_assemblies.organellar_genome,
        MAIN_WORKFLOW_VALIDATE_TAXID.out.versions,
        params.include,
        params.exclude,
        fcs_gx_database_path,
        ch_grabbed_reads_path
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
            MAIN_WORKFLOW_VALIDATE_TAXID.out.versions,
            organellar_include,
            organellar_exclude,
            fcs_gx_database_path,
            ch_grabbed_reads_path
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
        params.hook_url
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
