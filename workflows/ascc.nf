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
include { VALIDATE_TAXID                                } from '../modules/local/validate_taxid'
include { FILTER_FASTA                                  } from '../modules/local/filter_fasta'
include { CREATE_BTK_DATASET                            } from '../modules/local/create_btk_dataset'
include { MERGE_BTK_DATASETS                            } from '../modules/local/merge_btk_datasets'
include { ASCC_MERGE_TABLES                             } from '../modules/local/ascc_merge_tables'
include { AUTOFILTER_AND_CHECK_ASSEMBLY                 } from '../modules/local/autofiltering'
include { SANGER_TOL_BTK                                } from '../modules/local/sanger_tol_btk'
include { GENERATE_SAMPLESHEET                          } from '../modules/local/generate_samplesheet'

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
    ch_versions             = Channel.empty()
    ch_out_merge            = Channel.empty()

    include_workflow_steps  = params.include ? params.include.split(",") : ""
    exclude_workflow_steps  = params.exclude ? params.exclude.split(",") : ""

    btk_busco_run_mode      = params.btk_busco_run_mode ? params.btk_busco_run_mode : "conditional"

    full_list               = ["kmers", "tiara", "coverage", "nt_blast", "nr_diamond", "uniprot_diamond", "kraken", "fcs-gx", "fcs-adaptor", "vecscreen", "btk_busco", "pacbio_barcodes", "organellar_blast", "autofilter_assembly", "ALL", ""]

    if (!full_list.containsAll(include_workflow_steps) && !full_list.containsAll(exclude_workflow_steps)) {
        exit 1, "There is an extra argument given on Command Line: \n Check contents of: $include_workflow_steps\nAnd $exclude_workflow_steps\nMaster list is: $full_list"
    }

    input_ch                = Channel.fromPath(params.input, checkIfExists: true)

    //
    // SUBWORKFLOW: DECODE YAML INTO PARAMETERS FOR PIPELINE
    //
    YAML_INPUT (
        input_ch
    )
    ch_versions             = ch_versions.mix(YAML_INPUT.out.versions)

    //
    // MODULE: ENSURE THAT THE TAXID FOR THE INPUT GENOME IS INDEED IN THE TAXDUMP
    //
    VALIDATE_TAXID(
        YAML_INPUT.out.taxid,
        YAML_INPUT.out.ncbi_taxonomy_path
    )
    ch_versions             = ch_versions.mix(VALIDATE_TAXID.out.versions)


    //
    // MODULE:
    //
    FILTER_FASTA(
        YAML_INPUT.out.reference_tuple,
    )
    ch_versions             = ch_versions.mix(FILTER_FASTA.out.versions)


    //
    // LOGIC: INJECT SLIDING WINDOW VALUES INTO REFERENCE
    //
    FILTER_FASTA.out.fasta
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
    // MODULE: CALCULATE GC CONTENT PER SCAFFOLD IN INPUT FASTA
    //
    GC_CONTENT (
        FILTER_FASTA.out.fasta
    )
    ch_versions             = ch_versions.mix(GC_CONTENT.out.versions)

    //
    // SUBWORKFLOW: GENERATE GENOME FILE
    //
    GENERATE_GENOME (
        FILTER_FASTA.out.fasta,
        YAML_INPUT.out.pacbio_barcodes
    )
    ch_versions             = ch_versions.mix(GENERATE_GENOME.out.versions)

    //
    // SUBWORKFLOW: GENERATE A REPORT ON LENGTHS OF N's IN THE INPUT GENOMe
    //
    TRAILINGNS_CHECK (
        FILTER_FASTA.out.fasta
    )
    ch_versions             = ch_versions.mix(TRAILINGNS_CHECK.out.versions)

    //
    // SUBWORKFLOW: COUNT KMERS, THEN REDUCE DIMENSIONS USING SELECTED METHODS
    //

    if ( include_workflow_steps.contains('kmers') || include_workflow_steps.contains('ALL')) {

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
        ch_versions         = ch_versions.mix(GET_KMERS_PROFILE.out.versions)
        ch_kmers            = GET_KMERS_PROFILE.out.combined_csv
    } else {
        ch_kmers            = []
    }

    //
    // SUBWORKFLOW: EXTRACT RESULTS HITS FROM TIARA
    //
    if ( include_workflow_steps.contains('tiara') || include_workflow_steps.contains('ALL')) {
        EXTRACT_TIARA_HITS (
            GENERATE_GENOME.out.reference_tuple
        )
        ch_versions         = ch_versions.mix(EXTRACT_TIARA_HITS.out.versions)
        ch_tiara            = EXTRACT_TIARA_HITS.out.ch_tiara.map{it[1]}
    } else {
        ch_tiara            = []
    }

    //
    // SUBWORKFLOW: EXTRACT RESULTS HITS FROM NT-BLAST
    //
    if ( include_workflow_steps.contains('nt_blast') || include_workflow_steps.contains('ALL') ) {
        //
        // NOTE: ch_nt_blast needs to be set in two places incase it
        //          fails during the run
        //

        ch_nt_blast         = []
        EXTRACT_NT_BLAST (
            modified_input,
            YAML_INPUT.out.nt_database,
            YAML_INPUT.out.ncbi_accessions,
            YAML_INPUT.out.ncbi_rankedlineage_path
        )
        ch_versions         = ch_versions.mix(EXTRACT_NT_BLAST.out.versions)
        ch_nt_blast         = EXTRACT_NT_BLAST.out.ch_blast_hits.map{it[1]}

    } else {
        ch_nt_blast         = []
    }

    if ( include_workflow_steps.contains('organellar_blast') || include_workflow_steps.contains('ALL') ) {
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
            GENERATE_GENOME.out.reference_tuple,
            YAML_INPUT.out.mito_var,
            mito_check.valid
        )
        ch_mito             = MITO_ORGANELLAR_BLAST.out.organelle_report.map{it[1]}
        ch_versions         = ch_versions.mix(MITO_ORGANELLAR_BLAST.out.versions)
    } else {
        ch_mito             = []
    }

    if ( include_workflow_steps.contains('organellar_blast') || include_workflow_steps.contains('ALL') ) {

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
            GENERATE_GENOME.out.reference_tuple,
            YAML_INPUT.out.plastid_var,
            plastid_check.valid
        )
        ch_chloro           = PLASTID_ORGANELLAR_BLAST.out.organelle_report.map{it[1]}
        ch_versions         = ch_versions.mix(PLASTID_ORGANELLAR_BLAST.out.versions)
    } else {
        ch_chloro           = []
    }

    //
    // SUBWORKFLOW: RUN FCS-ADAPTOR TO IDENTIDY ADAPTOR AND VECTORR CONTAMINATION
    //
    if ( include_workflow_steps.contains('fcs-adaptor') || include_workflow_steps.contains('ALL') ) {
        RUN_FCSADAPTOR (
            GENERATE_GENOME.out.reference_tuple
        )
        RUN_FCSADAPTOR.out.ch_euk
            .map{it[1]}
            .combine(
                RUN_FCSADAPTOR.out.ch_prok.map{it[1]}
            )
            .set{ ch_fcsadapt }
        ch_versions         = ch_versions.mix(RUN_FCSADAPTOR.out.versions)
    } else {
        ch_fcsadapt         = []
    }

    //
    // SUBWORKFLOW: RUN FCS-GX TO IDENTIFY CONTAMINATION IN THE ASSEMBLY
    //
    if ( include_workflow_steps.contains('fcs-gx') || include_workflow_steps.contains('ALL') ) {
        RUN_FCSGX (
            GENERATE_GENOME.out.reference_tuple,
            YAML_INPUT.out.fcs_gx_database_path,
            YAML_INPUT.out.taxid,
            YAML_INPUT.out.ncbi_rankedlineage_path
        )

        ch_fcsgx            = RUN_FCSGX.out.fcsgxresult.map{it[1]}
        ch_versions         = ch_versions.mix(RUN_FCSGX.out.versions)
    } else {
        ch_fcsgx            = []
    }

    //
    // SUBWORKFLOW: IDENTITY PACBIO BARCODES IN INPUT DATA
    //
    if ( include_workflow_steps.contains('pacbio_barcodes') || include_workflow_steps.contains('ALL') ) {
        PACBIO_BARCODE_CHECK (
            GENERATE_GENOME.out.reference_tuple,
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
                ch_barcode // Not in use
            }

        ch_versions         = ch_versions.mix(PACBIO_BARCODE_CHECK.out.versions)
    } else {
        ch_barcode          = []
    }

    //
    // SUBWORKFLOW: CALCULATE AVERAGE READ COVERAGE
    //
    if ( include_workflow_steps.contains('coverage') || include_workflow_steps.contains('btk_busco') || include_workflow_steps.contains('ALL') ) {
        RUN_READ_COVERAGE (
            GENERATE_GENOME.out.reference_tuple,
            YAML_INPUT.out.assembly_path,
            YAML_INPUT.out.pacbio_tuple,
            YAML_INPUT.out.reads_type
        )
        ch_coverage         = RUN_READ_COVERAGE.out.tsv_ch.map{it[1]}
        ch_bam              = RUN_READ_COVERAGE.out.bam_ch.map{it[1]}
        ch_versions         = ch_versions.mix(RUN_READ_COVERAGE.out.versions)
    } else {
        ch_coverage         = []
        ch_bam              = []
    }

    //
    // SUBWORKFLOW: COLLECT SOFTWARE VERSIONS
    //
    if ( include_workflow_steps.contains('vecscreen') || include_workflow_steps.contains('ALL') ) {
        RUN_VECSCREEN (
            GENERATE_GENOME.out.reference_tuple,
            YAML_INPUT.out.vecscreen_database_path
        )
        ch_vecscreen        = RUN_VECSCREEN.out.vecscreen_contam.map{it[1]}
        ch_versions         = ch_versions.mix(RUN_VECSCREEN.out.versions)
    } else {
        ch_vecscreen        = []
    }

    //
    // SUBWORKFLOW: Run the kraken classifier
    //
    if ( include_workflow_steps.contains('kraken') || include_workflow_steps.contains('ALL') ) {
        RUN_NT_KRAKEN(
            GENERATE_GENOME.out.reference_tuple,
            YAML_INPUT.out.nt_kraken_db_path,
            YAML_INPUT.out.ncbi_rankedlineage_path
        )
        ch_kraken1          = RUN_NT_KRAKEN.out.classified.map{it[1]}
        ch_kraken2          = RUN_NT_KRAKEN.out.report.map{it[1]}
        ch_kraken3          = RUN_NT_KRAKEN.out.lineage

        ch_versions         = ch_versions.mix(RUN_NT_KRAKEN.out.versions)
    } else {
        ch_kraken1          = []
        ch_kraken2          = []
        ch_kraken3          = []
    }

    //
    // SUBWORKFLOW: DIAMOND BLAST FOR INPUT ASSEMBLY
    //
    if ( include_workflow_steps.contains('nr_diamond') || include_workflow_steps.contains('ALL') ) {
        NUCLEOT_DIAMOND (
            modified_input,
            YAML_INPUT.out.diamond_nr_database_path
        )
        nr_full             = NUCLEOT_DIAMOND.out.reformed.map{it[1]}
        nr_hits             = NUCLEOT_DIAMOND.out.hits_file.map{it[1]}
        ch_versions         = ch_versions.mix(NUCLEOT_DIAMOND.out.versions)
    } else {
        nr_hits             = []
        nr_full             = []
    }

    //
    // SUBWORKFLOW: DIAMOND BLAST FOR INPUT ASSEMBLY
    //
    //qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms sphylums salltitles
    if ( include_workflow_steps.contains('uniprot_diamond') || include_workflow_steps.contains('ALL') ) {
        UNIPROT_DIAMOND (
            modified_input,
            YAML_INPUT.out.diamond_uniprot_database_path
        )
        un_full             = UNIPROT_DIAMOND.out.reformed.map{it[1]}
        un_hits             = UNIPROT_DIAMOND.out.hits_file.map{it[1]}
        ch_versions         = ch_versions.mix(UNIPROT_DIAMOND.out.versions)
    } else {
        un_hits             = []
        un_full             = []
    }

    ch_dot_genome           = GENERATE_GENOME.out.dot_genome.map{it[1]}

    CREATE_BTK_DATASET (
        GENERATE_GENOME.out.reference_tuple,
        ch_dot_genome,
        ch_kmers,
        ch_tiara,
        ch_nt_blast,
        ch_fcsgx,
        ch_bam,
        ch_coverage,
        ch_kraken1,
        ch_kraken2,
        ch_kraken3,
        nr_full,
        un_full,
        YAML_INPUT.out.ncbi_taxonomy_path.first()
    )
    ch_versions             = ch_versions.mix(CREATE_BTK_DATASET.out.versions)


    //
    // MODULE: AUTOFILTER ASSEMBLY BY TIARA AND FCSGX RESULTS
    //
    if ( include_workflow_steps.contains('tiara') && include_workflow_steps.contains('fcs-gx') && include_workflow_steps.contains("autofilter_assembly") || include_workflow_steps.contains('ALL') ) {
        AUTOFILTER_AND_CHECK_ASSEMBLY (
            GENERATE_GENOME.out.reference_tuple,
            EXTRACT_TIARA_HITS.out.ch_tiara,
            RUN_FCSGX.out.fcsgxresult,
            YAML_INPUT.out.ncbi_rankedlineage_path
        )
        ch_autofilt_assem   = AUTOFILTER_AND_CHECK_ASSEMBLY.out.decontaminated_assembly.map{it[1]}
        ch_autofilt_indicator = AUTOFILTER_AND_CHECK_ASSEMBLY.out.indicator_file

        AUTOFILTER_AND_CHECK_ASSEMBLY.out.alarm_file
            .map { file -> file.text.trim() }
            .branch { it ->
                run_btk: "ABNORMAL" ? it.contains("YES_ABNORMAL"): false
                dont_run: []
            }
            .set { btk_bool }


        ch_versions         = ch_versions.mix(AUTOFILTER_AND_CHECK_ASSEMBLY.out.versions)
    } else {
        ch_autofilt_assem   = []
        ch_autofilt_indicator = []
    }

    //
    // PIPELINE: PREPARE THE DATA FOR USE IN THE SANGER-TOL/BLOBTOOLKIT PIPELINE
    //              WE ARE USING THE PIPELINE HERE AS A MODULE THIS REQUIRES IT
    //              TO BE USED AS A AN INTERACTIVE JOB ON WHAT EVER EXECUTOR YOU ARE USING.
    //              This will also eventually check for the above run_btk boolean from autofilter
    if ( !exclude_workflow_steps.contains("btk_busco") && include_workflow_steps.contains('btk_busco') && btk_busco_run_mode == "conditional" && include_workflow_steps.contains("autofilter_assembly") && btk_bool.run_btk == "ABNORMAL" || !exclude_workflow_steps.contains("btk_busco") && include_workflow_steps.contains('ALL') || btk_busco_run_mode == "mandatory" && !exclude_workflow_steps.contains('btk_busco') && include_workflow_steps.contains('btk_busco') ) {
        GENERATE_GENOME.out.reference_tuple
            .combine(ch_bam)
            .map{ meta, ref, bam ->
                tuple(  [   id: meta.id ],
                        bam
                )
            }
            .set { new_bam }

        GENERATE_SAMPLESHEET (
            new_bam
        )
        ch_versions         = ch_versions.mix(GENERATE_SAMPLESHEET.out.versions)

        SANGER_TOL_BTK (
            GENERATE_GENOME.out.reference_tuple,
            new_bam,
            GENERATE_SAMPLESHEET.out.csv,
            YAML_INPUT.out.diamond_uniprot_database_path,
            YAML_INPUT.out.nt_database.map{it -> it[1]},
            YAML_INPUT.out.diamond_uniprot_database_path,
            [],
            YAML_INPUT.out.ncbi_taxonomy_path,
            YAML_INPUT.out.btk_yaml,
            YAML_INPUT.out.busco_lineages,
            YAML_INPUT.out.taxid,
            'GCA_0001'
        )
        ch_versions         = ch_versions.mix(SANGER_TOL_BTK.out.versions)


        MERGE_BTK_DATASETS (
            CREATE_BTK_DATASET.out.btk_datasets,
            SANGER_TOL_BTK.out.dataset
        )
        ch_versions         = ch_versions.mix(MERGE_BTK_DATASETS.out.versions)
        busco_merge_btk     = MERGE_BTK_DATASETS.out.busco_summary_tsv.map{it[1]}

    } else {
        busco_merge_btk     = []
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
        nr_hits,                                            // FROM -- NUCLEOT_DIAMOND.out.reformed.map{it[1]}
        un_hits,                                            // FROM -- UNIPROT_DIAMOND.out.reformed.map{it[1]}
        [],                                                 // <-- MARKER SCAN -- NOT IN PIPELINE YET
        [],                                                 // <-- CONTIGVIZ -- NOT IN PIPELINE YET
        CREATE_BTK_DATASET.out.create_summary.map{it[1]},   // FROM -- CREATE_BTK_DATASET
        busco_merge_btk,                                    // FROM -- MERGE_BTK_DATASETS.out.busco_summary_tsv
        ch_fcsgx                                            // FROM -- PARSE_FCSGX_RESULT.out.fcsgxresult.map{it[1]}
    )
    ch_versions             = ch_versions.mix(ASCC_MERGE_TABLES.out.versions)


    GenIndicator (
        ch_autofilt_indicator,
        ch_fcsgx,
        ch_fcsadapt,
        ch_tiara,
        ch_vecscreen,
        ch_barcode,
    )


    //
    // SUBWORKFLOW: Collates version data from prior subworflows
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    emit:
    software_ch             = CUSTOM_DUMPSOFTWAREVERSIONS.out.yml
    versions_ch             = CUSTOM_DUMPSOFTWAREVERSIONS.out.versions
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

process GenIndicator {
    label 'process_tiny'

    tag "Generating Phase 1 Indicator"
    executor 'local'

    input:
    val(a)
    val(b)
    val(c)
    val(d)
    val(e)
    val(f)

    output:
    path("decon_first_stage_done_indicator_file.txt")

    script:
    "touch decon_first_stage_done_indicator_file.txt"
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

    if (workflow.success) {
        // Generate a pipeline completion indicator file, for use in
        // Sanger-ToL automation
        def newFile = new File("${params.outdir}/pipeline_run_done_indicator_file.txt")
        newFile.createNewFile()
    }

}

workflow.onError = {
    def newFile = new File("${params.outdir}/pipeline_run_ERROR_indicator_file.txt")
    newFile.createNewFile()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
