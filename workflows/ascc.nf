/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowAscc.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { YAML_INPUT           } from '../subworkflows/local/yaml_input'
include { GENERATE_GENOME      } from '../subworkflows/local/generate_genome'
include { EXTRACT_TIARA_HITS   } from '../subworkflows/local/extract_tiara_hits'
include { EXTRACT_NT_BLAST     } from '../subworkflows/local/extract_nt_blast'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASCC {

    main:
    ch_versions = Channel.empty()

    input_ch = Channel.fromPath(params.input, checkIfExists: true)

    //
    // SUBWORKFLOW: DECODE YAML INTO PARAMETERS FOR PIPELINE
    //
    YAML_INPUT (
        input_ch
    )
    ch_versions = ch_versions.mix(YAML_INPUT.out.versions)

    //Channel.fromPath( YAML_INPUT.out.nt_database, checkIfExists=true )
//        .set { blast_db }

    //
    // SUBWORKFLOW: GENERATE GENOME FILE
    //
    GENERATE_GENOME (
        YAML_INPUT.out.assembly_title,
        YAML_INPUT.out.reference
    )
    ch_versions = ch_versions.mix(GENERATE_GENOME.out.versions)

    //
    // SUBWORKFLOW: EXTRACT RESULTS HITS FROM TIARA
    //
    EXTRACT_TIARA_HITS (
        GENERATE_GENOME.out.reference_tuple
    )
    ch_versions = ch_versions.mix(EXTRACT_TIARA_HITS.out.versions)

    //
    // LOGIC: INJECT SLIDING WINDOW VALUES INTO REFERENCE
    //
    GENERATE_GENOME.out.reference_tuple
        .combine ( YAML_INPUT.out.seqkit_sliding.toInteger() )
        .combine ( YAML_INPUT.out.seqkit_window.toInteger() )
        .map { meta, ref, sliding, window ->
            tuple([ id      : meta.id,
                    sliding : sliding,
                    window  : window
                ],
                file(ref)
            )}
        .set { modified_input }

    //
    // SUBWORKFLOW: EXTRACT RESULTS HITS FROM NT-BLAST
    //
    EXTRACT_NT_BLAST (
        modified_input,
        YAML_INPUT.out.nt_database,
        YAML_INPUT.out.ncbi_taxonomy_path,
        YAML_INPUT.out.ncbi_rankedlineage_path
    )
    ch_versions = ch_versions.mix(EXTRACT_NT_BLAST.out.versions)

    //
    // SUBWORKFLOW: COLLECT SOFTWARE VERSIONS
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    emit:
    software_ch = CUSTOM_DUMPSOFTWAREVERSIONS.out.yml
    versions_ch = CUSTOM_DUMPSOFTWAREVERSIONS.out.versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
