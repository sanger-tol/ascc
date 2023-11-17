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

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { YAML_INPUT                                    } from '../subworkflows/local/yaml_input'
include { GENERATE_GENOME                               } from '../subworkflows/local/generate_genome'
include { EXTRACT_TIARA_HITS                            } from '../subworkflows/local/extract_tiara_hits'
include { EXTRACT_NT_BLAST                              } from '../subworkflows/local/extract_nt_blast'
include { RUN_FCSADAPTOR                                } from '../subworkflows/local/run_fcsadaptor'
include { RUN_NT_KRAKEN                                 } from '../subworkflows/local/run_nt_kraken'
include { RUN_FCSGX                                     } from '../subworkflows/local/run_fcsgx'
include { PACBIO_BARCODE_CHECK                          } from '../subworkflows/local/pacbio_barcode_check'
include { RUN_READ_COVERAGE                             } from '../subworkflows/local/run_read_coverage'
include { RUN_VECSCREEN                                 } from '../subworkflows/local/run_vecscreen'
include { ORGANELLAR_BLAST as PLASTID_ORGANELLAR_BLAST  } from '../subworkflows/local/organellar_blast'
include { ORGANELLAR_BLAST as MITO_ORGANELLAR_BLAST     } from '../subworkflows/local/organellar_blast'

//
// MODULE: Local modules
//
include { GC_CONTENT                                    } from '../modules/local/gc_content'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS                   } from '../modules/nf-core/custom/dumpsoftwareversions/main'

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

    //
    // MODULE: CALCULATE GC CONTENT PER SCAFFOLD IN INPUT FASTA
    //
    /*
    GC_CONTENT (
        YAML_INPUT.out.reference_tuple
    )
    ch_versions = ch_versions.mix(GC_CONTENT.out.versions)
    */

    //
    // SUBWORKFLOW: GENERATE GENOME FILE
    //
    GENERATE_GENOME (
        YAML_INPUT.out.reference_tuple,
        YAML_INPUT.out.pacbio_barcodes
    )
    ch_versions = ch_versions.mix(GENERATE_GENOME.out.versions)

    //
    // SUBWORKFLOW: EXTRACT RESULTS HITS FROM TIARA
    //

    /*
    EXTRACT_TIARA_HITS (
        GENERATE_GENOME.out.reference_tuple
    )
    ch_versions = ch_versions.mix(EXTRACT_TIARA_HITS.out.versions)
    */

    //
    // LOGIC: INJECT SLIDING WINDOW VALUES INTO REFERENCE
    //
    YAML_INPUT.out.reference_tuple
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
        YAML_INPUT.out.ncbi_accessions,
        YAML_INPUT.out.ncbi_rankedlineage_path
    )
    ch_versions = ch_versions.mix(EXTRACT_NT_BLAST.out.versions)

    //
    // LOGIC: CHECK WHETHER THERE IS A MITO AND BRANCH
    //
    YAML_INPUT.out.mito_tuple
        .branch { meta, check ->
            valid:      check != "NO MITO"
            invalid:    check == "NO MITO"
        }
        .set { mito_check }


    //
    // SUBWORKFLOW: BLASTING FOR MITO ASSEMBLIES IN GENOME
    //
    MITO_ORGANELLAR_BLAST (
        YAML_INPUT.out.reference_tuple,
        YAML_INPUT.out.mito_var,
        mito_check.valid
    )
    ch_versions = ch_versions.mix(MITO_ORGANELLAR_BLAST.out.versions)

    //
    // LOGIC: CHECK WHETHER THERE IS A PLASTID AND BRANCH
    //
    YAML_INPUT.out.plastid_tuple
        .branch { meta, check ->
            valid:      check != "NO PLASTID"
            invalid:    check == "NO PLASTID"
        }
        .set { plastid_check }

    //
    // SUBWORKFLOW: BLASTING FOR PLASTID ASSEMBLIES IN GENOME
    //
    PLASTID_ORGANELLAR_BLAST (
        YAML_INPUT.out.reference_tuple,
        YAML_INPUT.out.plastid_var,
        plastid_check.valid
    )
    ch_versions = ch_versions.mix(PLASTID_ORGANELLAR_BLAST.out.versions)

    //
    // SUBWORKFLOW:
    //
    RUN_FCSADAPTOR (
        YAML_INPUT.out.reference_tuple
    )
    ch_versions = ch_versions.mix(RUN_FCSADAPTOR.out.versions)

    //
    // SUBWORKFLOW:
    //
    RUN_FCSGX (
        YAML_INPUT.out.reference_tuple,
        YAML_INPUT.out.fcs_gx_database_path,
        YAML_INPUT.out.taxid,
        YAML_INPUT.out.ncbi_rankedlineage_path
    )
    ch_versions = ch_versions.mix(RUN_FCSADAPTOR.out.versions)

    //
    // SUBWORKFLOW: IDENTITY PACBIO BARCODES IN INPUT DATA
    //
    PACBIO_BARCODE_CHECK (
        YAML_INPUT.out.reference_tuple,
        YAML_INPUT.out.pacbio_tuple,
        YAML_INPUT.out.pacbio_barcodes,
        YAML_INPUT.out.pacbio_multiplex_codes
    )
    ch_versions = ch_versions.mix(PACBIO_BARCODE_CHECK.out.versions)

    //
    // SUBWORKFLOW: CALCULATE AVERAGE READ COVERAGE
    //
    RUN_READ_COVERAGE (
        YAML_INPUT.out.reference_tuple,
        YAML_INPUT.out.assembly_path,
        YAML_INPUT.out.pacbio_tuple,
        YAML_INPUT.out.reads_type
    )
    ch_versions = ch_versions.mix(RUN_READ_COVERAGE.out.versions)

    //
    // SUBWORKFLOW: COLLECT SOFTWARE VERSIONS
    //
    RUN_VECSCREEN (
        GENERATE_GENOME.out.reference_tuple,
        YAML_INPUT.out.vecscreen_database_path
    )
    ch_versions = ch_versions.mix(RUN_VECSCREEN.out.versions)

    //
    // SUBWORKFLOW: Collates version data from prior subworflows
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
    // TreeValProject.summary(workflow, reference_tuple, summary_params, projectDir)

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
