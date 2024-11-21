/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CREATE_BTK_DATASET                            } from '../modules/local/create_btk_dataset'
include { MERGE_BTK_DATASETS                            } from '../modules/local/merge_btk_datasets'
include { ASCC_MERGE_TABLES                             } from '../modules/local/ascc_merge_tables'
include { AUTOFILTER_AND_CHECK_ASSEMBLY                 } from '../modules/local/autofiltering'
include { SANGER_TOL_BTK                                } from '../modules/local/sanger_tol_btk'
include { GENERATE_SAMPLESHEET                          } from '../modules/local/generate_samplesheet'

include { ESSENTIAL_JOBS                                } from '../subworkflows/local/essential_jobs'
include { GET_KMERS_PROFILE                             } from '../subworkflows/local/get_kmers_profile'
include { EXTRACT_TIARA_HITS                            } from '../subworkflows/local/extract_tiara_hits'
include { EXTRACT_NT_BLAST                              } from '../subworkflows/local/extract_nt_blast'
include { ORGANELLAR_BLAST as PLASTID_ORGANELLAR_BLAST  } from '../subworkflows/local/organellar_blast'
include { ORGANELLAR_BLAST as MITO_ORGANELLAR_BLAST     } from '../subworkflows/local/organellar_blast'
include { PACBIO_BARCODE_CHECK                          } from '../subworkflows/local/pacbio_barcode_check'
include { TRAILINGNS_CHECK                              } from '../subworkflows/local/trailingns_check'
include { RUN_READ_COVERAGE                             } from '../subworkflows/local/run_read_coverage'
include { RUN_VECSCREEN                                 } from '../subworkflows/local/run_vecscreen'
include { RUN_NT_KRAKEN                                 } from '../subworkflows/local/run_nt_kraken'
include { RUN_FCSGX                                     } from '../subworkflows/local/run_fcsgx'
include { RUN_FCSADAPTOR                                } from '../subworkflows/local/run_fcsadaptor'
include { RUN_DIAMOND as NR_DIAMOND                     } from '../subworkflows/local/run_diamond.nf'
include { RUN_DIAMOND as UP_DIAMOND                     } from '../subworkflows/local/run_diamond.nf'

include { paramsSummaryMap                              } from 'plugin/nf-validation'
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
    ch_samplesheet // channel: samplesheet read in from --input
    organellar_genomes
    validate_taxid_versions // Versions channel from main.nf
    include_steps
    exclude_steps

    main:
    ch_versions = Channel.empty()
    ch_versions = ch_versions.mix(validate_taxid_versions)

    //
    // LOGIC: CONTROL OF THE INCLUDE AND EXCLUDE FLAGS
    //      TODO: THESE SHOULD CREATE A SET OF INCLUDE - EXCLUDE
    //
    include_workflow_steps  = include_steps ? include_steps.split(",") : "ALL"
    exclude_workflow_steps  = exclude_steps ? exclude_steps.split(",") : "NONE"

    full_list               = ["kmers", "tiara", "coverage", "nt_blast", "nr_diamond", "uniprot_diamond", "kraken", "fcs-gx", "fcs-adaptor", "vecscreen", "btk_busco", "pacbio_barcodes", "organellar_blast", "autofilter_assembly", "ALL", "NONE"]

    if (!full_list.containsAll(include_workflow_steps) && !full_list.containsAll(exclude_workflow_steps)) {
        exit 1, "There is an extra argument given on Command Line: \n Check contents of: $include_workflow_steps\nAnd $exclude_workflow_steps\nMaster list is: $full_list"
    }

    println "GENOMIC RUN -- INCLUDE STEPS INC.: $include_workflow_steps"
    println "GENOMIC RUN -- EXCLUDE STEPS INC.: $exclude_workflow_steps"


    //
    // LOGIC: CREATE btk_busco_run_mode VALUE
    //
    btk_busco_run_mode = params.btk_busco_run_mode ?: "conditional"


    //
    // SUBWORKFLOW: RUNS FILTER_FASTA, GENERATE .GENOME, CALCS GC_CONTENT AND FINDS RUNS OF N's
    //
    ESSENTIAL_JOBS(
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(ESSENTIAL_JOBS.out.versions)


    if ( include_workflow_steps.contains('kmers') && !exclude_workflow_steps.contains("kmers") || include_workflow_steps.contains('ALL') && !exclude_workflow_steps.contains("kmers")) {

        //
        // LOGIC: CONVERT THE CHANNEL INTO AN EPOCH COUNT FOR THE GET_KMER_PROFILE
        //
        ESSENTIAL_JOBS.out.reference_tuple_from_GG
            .map { meta, file ->
                file.countFasta() * 3
            }
            .set {autoencoder_epochs_count}

        //
        // SUBWORKFLOW: COUNT KMERS, THEN REDUCE DIMENSIONS USING SELECTED METHODS
        //
        GET_KMERS_PROFILE (
            ESSENTIAL_JOBS.out.reference_tuple_from_GG,
            params.kmer_length,
            params.dimensionality_reduction_methods,
            autoencoder_epochs_count
        )
        ch_versions         = ch_versions.mix(GET_KMERS_PROFILE.out.versions)
        ch_kmers            = GET_KMERS_PROFILE.out.combined_csv
    } else {
        ch_kmers            = []
    }


    //
    // SUBWORKFLOW: EXTRACT RESULTS HITS FROM TIARA
    //
    if ( include_workflow_steps.contains('tiara') && !exclude_workflow_steps.contains("tiara") || include_workflow_steps.contains('ALL') && !exclude_workflow_steps.contains("tiara") ) {
        EXTRACT_TIARA_HITS (
            ESSENTIAL_JOBS.out.reference_tuple_from_GG
        )
        ch_versions         = ch_versions.mix(EXTRACT_TIARA_HITS.out.versions)
        ch_tiara            = EXTRACT_TIARA_HITS.out.ch_tiara.map{it[1]}
    } else {
        ch_tiara            = []
    }


    //
    // SUBWORKFLOW: EXTRACT RESULTS HITS FROM NT-BLAST
    //
    if ( include_workflow_steps.contains('nt_blast') && !exclude_workflow_steps.contains("nt_blast") || include_workflow_steps.contains('ALL') && !exclude_workflow_steps.contains("nt_blast")) {
        //
        // NOTE: ch_nt_blast needs to be set in two places incase it
        //          fails during the run
        //
        ch_nt_blast         = []
        ch_blast_lineage    = []

        EXTRACT_NT_BLAST (
            ESSENTIAL_JOBS.out.reference_tuple_from_GG,
            params.nt_database_path,
            params.ncbi_accession_ids_folder,
            params.ncbi_ranked_lineage_path
        )
        ch_versions         = ch_versions.mix(EXTRACT_NT_BLAST.out.versions)
        ch_nt_blast         = EXTRACT_NT_BLAST.out.ch_blast_hits.map{it[1]}
        ch_blast_lineage    = EXTRACT_NT_BLAST.out.ch_top_lineages.map{it[1]}

    } else {
        ch_nt_blast         = []
        ch_blast_lineage    = []
    }


    if ( include_workflow_steps.contains('organellar_blast') || include_workflow_steps.contains('ALL') ) {

        //
        // LOGIC: CHECK WHETHER THERE IS A MITO AND BRANCH
        //
        organellar_genomes
            .branch { meta, assembly ->
                mito:       meta.assembly_type == "MITO"
                plastid:    meta.assembly_type == "PLASTID"
                invalid:    true    // if value but not of the above conditions
            }
            .set { organellar_check }

        //
        // SUBWORKFLOW: BLASTING FOR MITO ASSEMBLIES IN GENOME
        //
        MITO_ORGANELLAR_BLAST (
            ESSENTIAL_JOBS.out.reference_tuple_from_GG,
            organellar_check.mito
        )

        ch_mito             = MITO_ORGANELLAR_BLAST.out.organelle_report.map{it[1]}
        ch_versions         = ch_versions.mix(MITO_ORGANELLAR_BLAST.out.versions)

        //
        // SUBWORKFLOW: BLASTING FOR PLASTID ASSEMBLIES IN GENOME
        //
        PLASTID_ORGANELLAR_BLAST (
            ESSENTIAL_JOBS.out.reference_tuple_from_GG,
            organellar_check.plastid
        )
        ch_chloro           = PLASTID_ORGANELLAR_BLAST.out.organelle_report.map{it[1]}
        ch_versions         = ch_versions.mix(PLASTID_ORGANELLAR_BLAST.out.versions)

    } else {
        ch_mito             = []
        ch_chloro           = []
    }


    //
    // SUBWORKFLOW: IDENTITY PACBIO BARCODES IN INPUT DATA
    //
    if ( include_workflow_steps.contains('pacbio_barcodes') || include_workflow_steps.contains('ALL') ) {
        PACBIO_BARCODE_CHECK (
            ESSENTIAL_JOBS.out.reference_tuple_from_GG,
            params.reads_path,
            params.reads_type,
            params.pacbio_barcode_file,
            params.pacbio_barcode_names
        )

        ch_versions         = ch_versions.mix(PACBIO_BARCODE_CHECK.out.versions)
    }


    //
    // SUBWORKFLOW: RUN FCS-ADAPTOR TO IDENTIDY ADAPTOR AND VECTORR CONTAMINATION
    //
    if ( include_workflow_steps.contains('fcs-adaptor') || include_workflow_steps.contains('ALL') ) {
        RUN_FCSADAPTOR (
            ESSENTIAL_JOBS.out.reference_tuple_from_GG
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
            ESSENTIAL_JOBS.out.reference_tuple_from_GG,
            params.fcs_gx_database_path,
            params.taxid,
            params.ncbi_ranked_lineage_path
        )

        ch_fcsgx            = RUN_FCSGX.out.fcsgxresult.map{it[1]}
        ch_versions         = ch_versions.mix(RUN_FCSGX.out.versions)
    } else {
        ch_fcsgx            = []
    }


    //
    // SUBWORKFLOW: CALCULATE AVERAGE READ COVERAGE
    //
    if ( include_workflow_steps.contains('coverage') || include_workflow_steps.contains('btk_busco') || include_workflow_steps.contains('ALL') ) {
        RUN_READ_COVERAGE (
            ESSENTIAL_JOBS.out.reference_tuple_from_GG,
            params.reads_path,
            params.reads_type,
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
            ESSENTIAL_JOBS.out.reference_tuple_from_GG,
            params.vecscreen_database_path
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
            ESSENTIAL_JOBS.out.reference_tuple_from_GG,
            params.nt_kraken_database_path,
            params.ncbi_ranked_lineage_path
        )
        ch_kraken1          = RUN_NT_KRAKEN.out.classified.map{it[1]}
        ch_kraken2          = RUN_NT_KRAKEN.out.report.map{it[1]}

        // TODO: Channel is not getting populated even though the it is includes.
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
    if ( include_workflow_steps.contains('nr_diamond') && !exclude_workflow_steps.contains("nr_diamond") || include_workflow_steps.contains('ALL') && !exclude_workflow_steps.contains("nr_diamond")) {
        NR_DIAMOND (
            ESSENTIAL_JOBS.out.reference_tuple_from_GG,
            params.diamond_nr_database_path
        )
        nr_full             = NR_DIAMOND.out.reformed.map{it[1]}
        nr_hits             = NR_DIAMOND.out.hits_file.map{it[1]}
        ch_versions         = ch_versions.mix(NR_DIAMOND.out.versions)
    } else {
        nr_hits             = []
        nr_full             = []
    }


    //
    // SUBWORKFLOW: DIAMOND BLAST FOR INPUT ASSEMBLY
    //  qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms sphylums salltitles
    if ( include_workflow_steps.contains('uniprot_diamond') || include_workflow_steps.contains('ALL') ) {
        UP_DIAMOND (
            ESSENTIAL_JOBS.out.reference_tuple_from_GG,
            params.diamond_uniprot_database_path
        )
        un_full             = UP_DIAMOND.out.reformed.map{it[1]}
        un_hits             = UP_DIAMOND.out.hits_file.map{it[1]}
        ch_versions         = ch_versions.mix(UP_DIAMOND.out.versions)
    } else {
        un_hits             = []
        un_full             = []
    }

    ch_dot_genome           = ESSENTIAL_JOBS.out.dot_genome.map{it[1]}


    //
    // MODULE: CREATE A BTK COMPATIBLE DATASET FOR NEW DATA
    //
    CREATE_BTK_DATASET (
        ESSENTIAL_JOBS.out.reference_tuple_from_GG,
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
        Channel.of(params.ncbi_taxonomy_path).first()
    )
    ch_versions             = ch_versions.mix(CREATE_BTK_DATASET.out.versions)


    //
    // MODULE: AUTOFILTER ASSEMBLY BY TIARA AND FCSGX RESULTS
    //
    if ( include_workflow_steps.contains('tiara') && include_workflow_steps.contains('fcs-gx') && include_workflow_steps.contains("autofilter_assembly") || include_workflow_steps.contains('ALL') ) {
        //
        // LOGIC: FILTER THE INPUT FOR THE AUTOFILTER STEP
        //          - We can't just combine on meta.id as some of the Channels have other data in there too
        //              so we just sanitise, and _then_ combine on 0, and _then_ add back in the taxid as we
        //              need that for this process. Thankfully taxid is a param so easy enough to add back in.
        //
        ESSENTIAL_JOBS.out.reference_tuple_from_GG
            .map{meta, file ->
                tuple([id: meta.id], file)
            }
            .set{ ref_tuple }

        EXTRACT_TIARA_HITS.out.ch_tiara
            .map{meta, file ->
                tuple([id: meta.id], file)
            }
            .set{ tiara_tuple }

        RUN_FCSGX.out.fcsgxresult
            .map{meta, file ->
                tuple([id: meta.id], file)
            }
            .set{ fcs_tuple }

        ref_tuple
            .combine(tiara_tuple, by: 0) // Essentially merge on meta, which the above standardises
            .combine(fcs_tuple, by: 0)
            .combine(Channel.of(params.ncbi_ranked_lineage_path))
            .set{ auto_filt_input }

        //
        // LOGIC: NOW MULTIMAP THE CHANNELS INTO CONSTITUENT CHANNELS SO THAT WE CAN RUN THE AUTOFILTER
        //
        auto_filt_input
            .combine(Channel.of(params.taxid))
            .multiMap{
                meta, ref, tiara, fcs, ncbi, thetaxid ->
                    reference: tuple([id: meta.id, taxid: thetaxid], ref)
                    tiara_file: tuple(meta, tiara)
                    fcs_file: tuple(meta, fcs)
                    ncbi_rank: ncbi
            }
            .set{ autofilter_input_formatted}

        //
        // MODULE: AUTOFILTER ASSEMBLY BY TIARA AND FCSGX RESULTS
        //
        AUTOFILTER_AND_CHECK_ASSEMBLY (
            autofilter_input_formatted.reference,
            autofilter_input_formatted.tiara_file,
            autofilter_input_formatted.fcs_file,
            autofilter_input_formatted.ncbi_rank
        )
        ch_autofilt_assem       = AUTOFILTER_AND_CHECK_ASSEMBLY.out.decontaminated_assembly.map{it[1]}
        ch_autofilt_indicator   = AUTOFILTER_AND_CHECK_ASSEMBLY.out.indicator_file

        //
        // LOGIC: BRANCH THE CHANNEL ON WHETHER OR NOT THERE IS ABNORMAL CONTAMINATION IN THE OUTPUT FILE.
        //
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


    if ( !exclude_workflow_steps.contains("btk_busco") && include_workflow_steps.contains('btk_busco') && btk_busco_run_mode == "conditional" && include_workflow_steps.contains("autofilter_assembly") && btk_bool.run_btk == "ABNORMAL" || !exclude_workflow_steps.contains("btk_busco") && include_workflow_steps.contains('ALL') || btk_busco_run_mode == "mandatory" && !exclude_workflow_steps.contains('btk_busco') && include_workflow_steps.contains('btk_busco') ) {

        //
        // MODULE: THIS MODULE FORMATS THE INPUT DATA IN A SPECIFIC CSV FORMAT FOR USE IN THE BTK PIPELINE
        //
        GENERATE_SAMPLESHEET (
            RUN_READ_COVERAGE.out.bam_ch
        )
        ch_versions         = ch_versions.mix(GENERATE_SAMPLESHEET.out.versions)

        //
        // PIPELINE: PREPARE THE DATA FOR USE IN THE SANGER-TOL/BLOBTOOLKIT PIPELINE
        //              WE ARE USING THE PIPELINE HERE AS A MODULE THIS REQUIRES IT
        //              TO BE USED AS A AN INTERACTIVE JOB ON WHAT EVER EXECUTOR YOU ARE USING.
        //              This will also eventually check for the above run_btk boolean from autofilter
        SANGER_TOL_BTK (
            ESSENTIAL_JOBS.out.reference_tuple_from_GG,
            GENERATE_SAMPLESHEET.out.csv,
            params.diamond_uniprot_database_path,
            params.nt_database_path,
            params.diamond_uniprot_database_path,
            [],
            params.ncbi_taxonomy_path,
            params.btk_yaml,
            params.busco_lineages,
            params.taxid,
            'GCA_0001'
        )
        ch_versions         = ch_versions.mix(SANGER_TOL_BTK.out.versions)

        //
        // MODULE: MERGE THE TWO BTK FORMATTED DATASETS INTO ONE DATASET FOR EASIER USE
        //
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
        ESSENTIAL_JOBS.out.gc_content_txt,                  // FROM -- GC_COVERAGE.out.tsv
        ch_coverage,                                        // FROM -- RUN_COVERAGE.out.tsv.map{it[1]}
        ch_tiara,                                           // FROM -- TIARA_TIARA.out.classifications.map{it[1]}
        [],                                                 // <-- BACTERIAL KRAKEN -- NOT IN PIPELINE YET
        ch_kraken3,                                         // FROM -- RUN_NT_KRAKEN.out.lineage.map{it[1]}
        ch_blast_lineage,                                        // FROM -- EXTRACT_NT_BLAST.out.ch_blast_hits.map{it[1]}
        ch_kmers,                                           // FROM -- GET_KMERS_PROFILE.out.combined_csv
        nr_hits,                                            // FROM -- NR_DIAMOND.out.reformed.map{it[1]}
        un_hits,                                            // FROM -- UP_DIAMOND.out.reformed.map{it[1]}
        [],                                                 // <-- MARKER SCAN -- NOT IN PIPELINE YET
        [],                                                 // <-- CONTIGVIZ -- NOT IN PIPELINE YET
        CREATE_BTK_DATASET.out.create_summary.map{it[1]},   // FROM -- CREATE_BTK_DATASET
        busco_merge_btk,                                    // FROM -- MERGE_BTK_DATASETS.out.busco_summary_tsv
        ch_fcsgx                                            // FROM -- PARSE_FCSGX_RESULT.out.fcsgxresult.map{it[1]}
    )
    ch_versions             = ch_versions.mix(ASCC_MERGE_TABLES.out.versions)


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
