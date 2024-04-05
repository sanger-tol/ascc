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
include { GET_KMERS_PROFILE                             } from '../subworkflows/local/get_kmers_profile'
include { RUN_READ_COVERAGE                             } from '../subworkflows/local/run_read_coverage'
include { RUN_VECSCREEN                                 } from '../subworkflows/local/run_vecscreen'
include { ORGANELLAR_BLAST as PLASTID_ORGANELLAR_BLAST  } from '../subworkflows/local/organellar_blast'
include { ORGANELLAR_BLAST as MITO_ORGANELLAR_BLAST     } from '../subworkflows/local/organellar_blast'
include { RUN_DIAMOND as NUCLEOT_DIAMOND                } from '../subworkflows/local/run_diamond.nf'
include { RUN_DIAMOND as UNIPROT_DIAMOND                } from '../subworkflows/local/run_diamond.nf'

//
// MODULE: Local modules
//
include { GC_CONTENT                                    } from '../modules/local/gc_content'
include { CREATE_BTK_DATASET                            } from '../modules/local/create_btk_dataset'


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
    ch_versions     = Channel.empty()
    ch_out_merge    = Channel.empty()

    workflow_steps  = params.steps.split(",")

    input_ch        = Channel.fromPath(params.input, checkIfExists: true)

    //
    // SUBWORKFLOW: DECODE YAML INTO PARAMETERS FOR PIPELINE
    //
    YAML_INPUT (
        input_ch
    )
    ch_versions     = ch_versions.mix(YAML_INPUT.out.versions)

    //
    // MODULE: CALCULATE GC CONTENT PER SCAFFOLD IN INPUT FASTA
    //
    GC_CONTENT (
        YAML_INPUT.out.reference_tuple
    )
    ch_versions     = ch_versions.mix(GC_CONTENT.out.versions)

    //
    // SUBWORKFLOW: GENERATE GENOME FILE
    //
    GENERATE_GENOME (
        YAML_INPUT.out.reference_tuple,
        YAML_INPUT.out.pacbio_barcodes
    )
    ch_versions     = ch_versions.mix(GENERATE_GENOME.out.versions)

    //
    // SUBWORKFLOW: COUNT KMERS, THEN REDUCE DIMENSIONS USING SELECTED METHODS
    //

    if ( workflow_steps.contains('kmers') || workflow_steps.contains('ALL')) {

        GENERATE_GENOME.out.reference_tuple
            .map { meta, file ->
                tuple (
                    meta,
                    file,
                    file.countFasta() * 3
                )
            }
            .set {autoencoder_epochs_count}

        GET_KMERS_PROFILE (
            GENERATE_GENOME.out.reference_tuple,
            YAML_INPUT.out.kmer_len,
            YAML_INPUT.out.dimensionality_reduction_methods,
            YAML_INPUT.out.n_neighbours,
            autoencoder_epochs_count.map{it -> it[2]}
        )
        ch_versions     = ch_versions.mix(GET_KMERS_PROFILE.out.versions)
        ch_kmers        = GET_KMERS_PROFILE.out.combined_csv
    } else {
        ch_kmers        = []
    }

    //
    // SUBWORKFLOW: EXTRACT RESULTS HITS FROM TIARA
    //
    if ( workflow_steps.contains('tiara') || workflow_steps.contains('ALL')) {
        EXTRACT_TIARA_HITS (
            GENERATE_GENOME.out.reference_tuple
        )
        ch_versions     = ch_versions.mix(EXTRACT_TIARA_HITS.out.versions)
        ch_tiara        = EXTRACT_TIARA_HITS.out.ch_tiara.map{it[1]}
    } else {
        ch_tiara        = []
    }

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
    if ( workflow_steps.contains('nt_blast') || workflow_steps.contains('ALL') ) {
        EXTRACT_NT_BLAST (
            modified_input,
            YAML_INPUT.out.nt_database,
            YAML_INPUT.out.ncbi_accessions,
            YAML_INPUT.out.ncbi_rankedlineage_path
        )
        ch_versions     = ch_versions.mix(EXTRACT_NT_BLAST.out.versions)
        ch_nt_blast     = EXTRACT_NT_BLAST.out.ch_blast_hits.map{it[1]}
    } else {
        ch_nt_blast     = []
    }

    if ( workflow_steps.contains('mito') || workflow_steps.contains('ALL') ) {
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
        ch_mito         = MITO_ORGANELLAR_BLAST.out.organelle_report.map{it[1]}
        ch_versions     = ch_versions.mix(MITO_ORGANELLAR_BLAST.out.versions)
    } else {
        ch_mito         = []
    }

    if ( workflow_steps.contains('chloro') || workflow_steps.contains('ALL') ) {

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
        ch_chloro       = PLASTID_ORGANELLAR_BLAST.out.organelle_report.map{it[1]}
        ch_versions     = ch_versions.mix(PLASTID_ORGANELLAR_BLAST.out.versions)
    } else {
        ch_chloro       = []
    }


    //
    // SUBWORKFLOW:
    //
    if ( workflow_steps.contains('fcs_adapt') || workflow_steps.contains('ALL') ) {
        RUN_FCSADAPTOR (
            YAML_INPUT.out.reference_tuple
        )
        RUN_FCSADAPTOR.out.ch_euk
            .map{it[1]}
            .combine(
                RUN_FCSADAPTOR.out.ch_prok.map{it[1]}
            )
            .set{ ch_fcsadapt }
        ch_fcsadapt
        ch_versions     = ch_versions.mix(RUN_FCSADAPTOR.out.versions)
    } else {
        ch_fcsadapt     = []
    }

    //
    // SUBWORKFLOW:
    //
    if ( workflow_steps.contains('fcsgx') || workflow_steps.contains('ALL') ) {
        RUN_FCSGX (
            YAML_INPUT.out.reference_tuple,
            YAML_INPUT.out.fcs_gx_database_path,
            YAML_INPUT.out.taxid,
            YAML_INPUT.out.ncbi_rankedlineage_path
        )
        ch_fcsgx        = RUN_FCSGX.out.fcsgxresult.map{it[1]}
        ch_versions     = ch_versions.mix(RUN_FCSADAPTOR.out.versions)
    } else {
        ch_fcsgx        = []
    }

    //
    // SUBWORKFLOW: IDENTITY PACBIO BARCODES IN INPUT DATA
    //
    if ( workflow_steps.contains('barcodes') || workflow_steps.contains('ALL') ) {
        PACBIO_BARCODE_CHECK (
            YAML_INPUT.out.reference_tuple,
            YAML_INPUT.out.pacbio_tuple,
            YAML_INPUT.out.pacbio_barcodes,
            YAML_INPUT.out.pacbio_multiplex_codes
        )

        PACBIO_BARCODE_CHECK.out.filtered
            .map{
                it[1]
            }
            .collect()
            .set {
                ch_barcode
            }

        ch_versions     = ch_versions.mix(PACBIO_BARCODE_CHECK.out.versions)
    } else {
        ch_barcode      = []
    }

    //
    // SUBWORKFLOW: CALCULATE AVERAGE READ COVERAGE
    //
    if ( workflow_steps.contains('coverage') || workflow_steps.contains('ALL') ) {
        RUN_READ_COVERAGE (
            YAML_INPUT.out.reference_tuple,
            YAML_INPUT.out.assembly_path,
            YAML_INPUT.out.pacbio_tuple,
            YAML_INPUT.out.reads_type
        )
        ch_coverage     = RUN_READ_COVERAGE.out.tsv_ch.map{it[1]}
        ch_bam          = RUN_READ_COVERAGE.out.bam_ch.map{it[1]}
        ch_versions     = ch_versions.mix(RUN_READ_COVERAGE.out.versions)
    } else {
        ch_coverage     = []
        ch_bam          = []
    }

    //
    // SUBWORKFLOW: COLLECT SOFTWARE VERSIONS
    //
    if ( workflow_steps.contains('vecscreen') || workflow_steps.contains('ALL') ) {
        RUN_VECSCREEN (
            GENERATE_GENOME.out.reference_tuple,
            YAML_INPUT.out.vecscreen_database_path
        )
        ch_vecscreen    = RUN_VECSCREEN.out.vecscreen_contam.map{it[1]}
        ch_versions     = ch_versions.mix(RUN_VECSCREEN.out.versions)
    } else {
        ch_vecscreen    = []
    }

    //
    // SUBWORKFLOW: Run the kraken classifier
    //
    if ( workflow_steps.contains('kraken') || workflow_steps.contains('ALL') ) {
        RUN_NT_KRAKEN(
            GENERATE_GENOME.out.reference_tuple,
            YAML_INPUT.out.nt_kraken_db_path,
            YAML_INPUT.out.ncbi_rankedlineage_path
        )
        ch_kraken1      = RUN_NT_KRAKEN.out.classifiedreformed
        ch_kraken2      = RUN_NT_KRAKEN.out.report.map{it[1]}
        ch_kraken3      = RUN_NT_KRAKEN.out.lineage

        ch_versions     = ch_versions.mix(RUN_NT_KRAKEN.out.versions)
    } else {
        ch_kraken1      = []
        ch_kraken2      = []
        ch_kraken3      = []
    }

    //
    // SUBWORKFLOW: DIAMOND BLAST FOR INPUT ASSEMBLY
    //
    if ( workflow_steps.contains('nt_diamond') || workflow_steps.contains('ALL') ) {
        NUCLEOT_DIAMOND (
            modified_input,
            YAML_INPUT.out.diamond_nr_database_path
        )
        nt_full         = NUCLEOT_DIAMOND.out.reformed.map{it[1]}
        ch_versions     = ch_versions.mix(RUN_DIAMOND.out.versions)
    } else {
        nt_full         = []
    }

    //
    // SUBWORKFLOW: DIAMOND BLAST FOR INPUT ASSEMBLY
    //
    if ( workflow_steps.contains('uniprot_diamond') || workflow_steps.contains('ALL') ) {
        UNIPROT_DIAMOND (
            modified_input,
            YAML_INPUT.out.diamond_uniprot_database_path
        )
        un_full         = UNIPROT_DIAMOND.out.reformed.map{it[1]}
        ch_versions     = ch_versions.mix(RUN_DIAMOND.out.versions)
    } else {
        un_full         = []
    }

    // mix the outputs of the outpuutting process so that we can
    // insert them into the one process to create the btk and the merged report
    // much like the versions channel

    CREATE_BTK_DATASET (
        GENERATE_GENOME.out.reference_tuple,
        GENERATE_GENOME.out.dot_genome.map{it[1]},
        ch_kmers,
        ch_tiara,
        ch_nt_blast,
        ch_fcsgx,
        ch_bam,
        ch_coverage,
        ch_kraken1,
        ch_kraken2,
        ch_kraken3,
        nt_full,
        un_full,
        YAML_INPUT.out.ncbi_taxonomy_path,

    )

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
