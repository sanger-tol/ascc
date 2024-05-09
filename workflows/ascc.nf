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
include { TRAILINGNS_CHECK                              } from '../subworkflows/local/trailingns_check'

//
// MODULE: Local modules
//
include { GC_CONTENT                                    } from '../modules/local/gc_content'

include { CREATE_BTK_DATASET                            } from '../modules/local/create_btk_dataset'
include { MERGE_BTK_DATASETS                            } from '../modules/local/merge_btk_datasets'
include { ASCC_MERGE_TABLES                             } from '../modules/local/ascc_merge_tables'
include { AUTOFILTER_AND_CHECK_ASSEMBLY                 } from '../modules/local/autofiltering'

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
    // SUBWORKFLOW: GENERATE A REPORT ON LENGTHS OF N's IN THE INPUT GENOMe
    //
    TRAILINGNS_CHECK (
        YAML_INPUT.out.reference_tuple
    )
    ch_versions = ch_versions.mix(TRAILINGNS_CHECK.out.versions)

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
        ch_kraken1      = RUN_NT_KRAKEN.out.classified.map{it[1]}
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
        nt_hits         = NUCLEOT_DIAMOND.out.hits_file.map{it[1]}
        ch_versions     = ch_versions.mix(NUCLEOT_DIAMOND.out.versions)
    } else {
        nt_hits         = []
        nt_full         = []
    }

    //
    // SUBWORKFLOW: DIAMOND BLAST FOR INPUT ASSEMBLY
    //
    //qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms sphylums salltitles
    if ( workflow_steps.contains('uniprot_diamond') || workflow_steps.contains('ALL') ) {
        UNIPROT_DIAMOND (
            modified_input,
            YAML_INPUT.out.diamond_uniprot_database_path
        )
        un_full         = UNIPROT_DIAMOND.out.reformed.map{it[1]}
        un_hits         = UNIPROT_DIAMOND.out.hits_file.map{it[1]}
        ch_versions     = ch_versions.mix(UNIPROT_DIAMOND.out.versions)
    } else {
        un_hits         = []
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
        nt_hits,
        un_hits,
        YAML_INPUT.out.ncbi_taxonomy_path,

    )
    ch_versions                 = ch_versions.mix(CREATE_BTK_DATASET.out.versions)


    //SANGER_TOL_BTK.out.btk_datasets = []
    //SANGER_TOL_BTK.out.summary = []


    //
    // MODULE: AUTOFILTER ASSEMBLY BY TIARA AND FCSGX RESULTS
    //
    if ( ( workflow_steps.contains('tiara') && workflow_steps.contains('fcsgx') && workflow_steps.contains("autofilter") ) || workflow_steps.contains('ALL') ) {
        AUTOFILTER_AND_CHECK_ASSEMBLY (
            YAML_INPUT.out.reference_tuple,
            EXTRACT_TIARA_HITS.out.ch_tiara,
            RUN_FCSGX.out.fcsgxresult
        )
        ch_autofiltered_assembly = AUTOFILTER_AND_CHECK_ASSEMBLY.out.decontaminated_assembly.map{it[1]}
        ch_versions              = ch_versions.mix(AUTOFILTER_AND_CHECK_ASSEMBLY.out.versions)
    } else {
        ch_autofiltered_assembly = []
    }

    //
    // LOGIC: SCAN FILE FOR PRESENCE OF ABNORMAL CONTAMINATION
    //          IF FOUND THEN WE WANT TO RUN BTK
    //
    ch_autofiltered_assembly
        .branch{
            btk_run: { if ch_autofiltered_assembly.getText().contains("YES_ABNORMAL_CONTAMINATION") ? "PASS" : [] }
            skip: []
        }
        .set { abnormal_flag }


    //
    // PIPELINE: PREPARE THE DATA FOR USE IN THE SANGER-TOL/BLOBTOOLKIT PIPELINE
    //              WE ARE USING THE PIPELINE HERE AS A MODULE THIS REQUIRES IT
    //              TO BE USED AS A AN INTERACTIVE JOB ON WHAT EVER EXECUTOR YOU ARE USING.
    //
    if ( ( workflow_steps.contains('busco_btk') && workflow_steps.contains("autofilter") && abnormal_flag ) || workflow_steps.contains('ALL') || workflow_steps.contains("force_btk") ) {

        GENERATE_SAMPLESHEET (
            YAML_INPUT.out.reference_tuple,
            YAML_INPUT.out.pacbio_tuple
        )
        ch_versions              = ch_versions.mix(GENERATE_SAMPLESHEET.out.versions)


        SANGER_TOL_BTK (
            YAML_INPUT.out.reference_tuple,
            GENERATE_SAMPLESHEET.out.csv,
            YAML_INPUT.out.blastp,
            YAML_INPUT.out.blastn,
            YAML_INPUT.out.blastx,
            [],
            YAML_INPUT.out.tax_dump,
            YAML_INPUT.out.taxon,
            'GCA_0001'
        )
        //ch_versions              = ch_versions.mix(SANGER_TOL_BTK.out.versions)


        MERGE_BTK_DATASETS (
            CREATE_BTK_DATASET.out.btk_datasets,
            [[],[]],     //SANGER_TOL_BTK.out.btk_datasets = []
            [[],[]]      //SANGER_TOL_BTK.out.summary = []
        )
        ch_versions              = ch_versions.mix(MERGE_BTK_DATASETS.out.versions)

    }


    //
    // SUBWORKFLOW: MERGES DATA THAT IS NOT USED IN THE CREATION OF THE BTK_DATASETS FOLDER
    //
    ASCC_MERGE_TABLES (
        GC_CONTENT.out.txt,                                 // FROM -- GC_COVERAGE.out.tsv
        ch_coverage,                                        // FROM -- RUN_COVERAGE.out.tsv.map{it[1]}
        ch_tiara,                                           // FROM -- TIARA_TIARA.out.classifications.map{it[1]}
        [],                                                 // <-- BACTERIAL KRAKEN -- NOT IN PIPELINE YET
        ch_kraken3,                                         // FROM -- RUN_NT_KRAKEN.out.lineage.map{it[1]}
        ch_nt_blast,                                        // FROM -- EXTRACT_NT_BLAST.out.ch_blast_hits.map{it[1]}
        ch_kmers,                                           // FROM -- GET_KMERS_PROFILE.out.combined_csv
        nt_hits,                                            // FROM -- NUCLEOT_DIAMOND.out.reformed.map{it[1]}
        un_hits,                                            // FROM -- UNIPROT_DIAMOND.out.reformed.map{it[1]}
        [],                                                 // <-- MARKER SCAN -- NOT IN PIPELINE YET
        [],                                                 // <-- CONTIGVIZ -- NOT IN PIPELINE YET
        CREATE_BTK_DATASET.out.create_summary.map{it[1]},
        [],                                                 // <-- BUSCO_BTK -- NOT IN PIPELINE YET
        ch_fcsgx                                            // FROM -- PARSE_FCSGX_RESULT.out.fcsgxresult.map{it[1]}
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

process GrabFiles {
    label 'process_tiny'

    tag "${meta.id}"
    executor 'local'

    input:
    tuple val(meta), path("in")

    output:
    tuple val(meta), path("in/*.{fa,fasta}.{gz}")

    "true"
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
