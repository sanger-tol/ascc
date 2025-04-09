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
include { NEXTFLOW_RUN as SANGER_TOL_BTK_CASCADE        } from '../modules/local/run/main'


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
include { GENERATE_HTML_REPORT_WORKFLOW                 } from '../subworkflows/local/generate_html_report'

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
    ch_samplesheet          // channel: samplesheet read in from --input
    organellar_genomes      // channel: tuple(meta, reference)
    validate_taxid_versions // Versions channel from main.nf
    include_steps           // params.include_steps
    exclude_steps           // params.exclude_steps
    fcs_db                  // [path(path)]
    reads

    main:
    ch_versions = Channel.empty()
    ch_versions = ch_versions.mix(validate_taxid_versions)

    //
    // LOGIC: CONTROL OF THE INCLUDE AND EXCLUDE FLAGS
    //      TODO: THESE SHOULD CREATE A SET OF INCLUDE - EXCLUDE
    //
    include_workflow_steps  = include_steps ? include_steps.split(",") : "ALL"
    exclude_workflow_steps  = exclude_steps ? exclude_steps.split(",") : "NONE"

    full_list               = [
        "essentials", "kmers", "tiara", "coverage", "nt_blast", "nr_diamond",
        "uniprot_diamond", "kraken", "fcs-gx", "fcs-adaptor", "vecscreen", "btk_busco",
        "pacbio_barcodes", "organellar_blast", "autofilter_assembly", "create_btk_dataset",
        "merge", "html_report", "ALL", "NONE"
    ]

    if (!full_list.containsAll(include_workflow_steps) && !full_list.containsAll(exclude_workflow_steps)) {
        exit 1, "There is an extra argument given on Command Line: \n Check contents of: $include_workflow_steps\nAnd $exclude_workflow_steps\nMaster list is: $full_list"
    }

    log.info "GENOMIC RUN -- INCLUDE STEPS INC.: $include_workflow_steps"
    log.info "GENOMIC RUN -- EXCLUDE STEPS INC.: $exclude_workflow_steps"


    //reads = CollectReads(reads_list)

    //
    // LOGIC: CREATE btk_busco_run_mode VALUE
    //
    btk_busco_run_mode = params.btk_busco_run_mode ?: "conditional"


    //
    // LOGIC: PRETTY NOTIFICATION OF FILES AT STAGE
    //
    ch_samplesheet
        .map { meta, sample ->
            log.info "GENOMIC WORKFLOW:\n\t-- $meta\n\t-- $sample"
        }


    //
    // SUBWORKFLOW: RUNS FILTER_FASTA, GENERATE .GENOME, CALCS GC_CONTENT AND FINDS RUNS OF N's
    //                  THIS SHOULD NOT RUN ONLY WHEN SPECIFICALLY REQUESTED
    //
    if ( !exclude_workflow_steps.contains("essentials")) {

        ESSENTIAL_JOBS(
            ch_samplesheet
        )
        ch_versions = ch_versions.mix(ESSENTIAL_JOBS.out.versions)

        reference_tuple_from_GG = ESSENTIAL_JOBS.out.reference_tuple_from_GG
        ej_dot_genome           = ESSENTIAL_JOBS.out.dot_genome
        ej_gc_coverage          = ESSENTIAL_JOBS.out.gc_content_txt
        reference_tuple_w_seqkt = ESSENTIAL_JOBS.out.reference_with_seqkit


    } else {
        log.warn("MAKE SURE YOU ARE AWARE YOU ARE SKIPPING ESSENTIAL JOBS, THIS INCLUDES BREAKING SCAFFOLDS OVER 1.9GB, FILTERING N\'s AND GC CONTENT REPORT (THIS WILL BREAK OTHER PROCESSES AND SHOULD ONLY BE RUN WITH `--include essentials`)")

        reference_tuple_from_GG = ch_samplesheet // This is the reference genome input channel
    }


    if ( (include_workflow_steps.contains('kmers') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("kmers")
    ) {
        //
        // LOGIC: CONVERT THE CHANNEL I AN EPOCH COUNT FOR THE GET_KMER_PROFILE
        //
        reference_tuple_from_GG
            .map { meta, file ->
                file.countFasta() * 3
            }
            .set {autoencoder_epochs_count}

        //
        // SUBWORKFLOW: COUNT KMERS, THEN REDUCE DIMENSIONS USING SELECTED METHODS
        //
        GET_KMERS_PROFILE (
            reference_tuple_from_GG,
            params.kmer_length,
            params.dimensionality_reduction_methods,
            autoencoder_epochs_count
        )
        ch_versions         = ch_versions.mix(GET_KMERS_PROFILE.out.versions)

        ch_kmers            = GET_KMERS_PROFILE.out.combined_csv
                                .map { it ->
                                    [[id: it[0].id, process: "KMERS"], it[1]]
                                }
                                .ifEmpty { [[],[]] }
    } else {
        ch_kmers            = Channel.of( [[],[]] )
    }


    //
    // SUBWORKFLOW: EXTRACT RESULTS HITS FROM TIARA
    //
    if ( (include_workflow_steps.contains('tiara') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("tiara")
    ) {
        EXTRACT_TIARA_HITS (
            reference_tuple_from_GG
        )
        ch_versions         = ch_versions.mix(EXTRACT_TIARA_HITS.out.versions)
        ch_tiara            = EXTRACT_TIARA_HITS.out.ch_tiara
                                .map { it ->
                                    [[id: it[0].id, process: "TIARA"], it[1]]
                                }
                                .ifEmpty { [[],[]] }
    } else {
        ch_tiara            = Channel.of( [[],[]] )
    }


    //
    // SUBWORKFLOW: EXTRACT RESULTS HITS FROM NT-BLAST
    //
    if ( (include_workflow_steps.contains('nt_blast') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("nt_blast")
    ) {
        //
        // NOTE: ch_nt_blast needs to be set in two places incase it
        //          fails during the run
        //
        EXTRACT_NT_BLAST (
            reference_tuple_from_GG,
            params.nt_database_path,
            params.ncbi_ranked_lineage_path
        )
        ch_versions         = ch_versions.mix(EXTRACT_NT_BLAST.out.versions)
        ch_nt_blast         = EXTRACT_NT_BLAST.out.ch_blast_hits
                                .map { it ->
                                    [[id: it[0].id, process: "NT-BLAST"], it[1]]
                                }
                                .ifEmpty { [[],[]] }

        ch_blast_lineage    = EXTRACT_NT_BLAST.out.ch_top_lineages
                                .map { it ->
                                    [[id: it[0].id, process: "NT-BLAST-LINEAGE"], it[1]]
                                }
                                .ifEmpty { [[],[]] }
                                
        ch_btk_format       = EXTRACT_NT_BLAST.out.ch_btk_format
                                .map { it ->
                                    [[id: it[0].id, process: "NT-BLAST-BTK"], it[1]]
                                }
                                .ifEmpty { [[],[]] }

    } else {
        ch_nt_blast         = Channel.of( [[],[]] )
        ch_blast_lineage    = Channel.of( [[],[]] )
        ch_btk_format       = Channel.of( [[],[]] )
    }


    //
    // SUBWORKFLOW: DIAMOND BLAST FOR INPUT ASSEMBLY
    //
    if ( (include_workflow_steps.contains('nr_diamond') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("nr_diamond")
    ) {
        NR_DIAMOND (
            reference_tuple_from_GG,
            params.diamond_nr_database_path
        )
        ch_versions         = ch_versions.mix(NR_DIAMOND.out.versions)
        nr_full             = NR_DIAMOND.out.reformed
                                .map { it ->
                                    [[id: it[0].id, process: "NR-FULL"], it[1]]
                                }
                                .ifEmpty { [[],[]] }

        nr_hits             = NR_DIAMOND.out.hits_file
                                .map { it ->
                                    [[id: it[0].id, process: "NR-HITS"], it[1]]
                                }
                                .ifEmpty { [[],[]] }

    } else {
        nr_full             = Channel.of( [[],[]] )
        nr_hits             = Channel.of( [[],[]] )
    }


    //
    // SUBWORKFLOW: DIAMOND BLAST FOR INPUT ASSEMBLY
    //  qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms sphylums salltitles
    if ( (include_workflow_steps.contains('uniprot_diamond') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("uniprot_diamond")
    ) {
        UP_DIAMOND (
            reference_tuple_from_GG,
            params.diamond_uniprot_database_path
        )
        ch_versions         = ch_versions.mix(UP_DIAMOND.out.versions)
        un_full             = UP_DIAMOND.out.reformed
                                .map { it ->
                                    [[id: it[0].id, process: "UN-FULL"], it[1]]
                                }
                                .ifEmpty { [[],[]] }

        un_hits             = UP_DIAMOND.out.hits_file
                                .map { it ->
                                    [[id: it[0].id, process: "UN-HITS"], it[1]]
                                }
                                .ifEmpty { [[],[]] }
    } else {
        un_full             = Channel.of( [[],[]] )
        un_hits             = Channel.of( [[],[]] )
    }


    if ( (include_workflow_steps.contains('organellar_blast') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("organellar_blast")
    ) {

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
            reference_tuple_from_GG,
            organellar_check.mito
        )
        ch_versions         = ch_versions.mix(MITO_ORGANELLAR_BLAST.out.versions)


        //
        // SUBWORKFLOW: BLASTING FOR PLASTID ASSEMBLIES IN GENOME
        //
        PLASTID_ORGANELLAR_BLAST (
            reference_tuple_from_GG,
            organellar_check.plastid
        )
        ch_versions         = ch_versions.mix(PLASTID_ORGANELLAR_BLAST.out.versions)

        ch_mito             = MITO_ORGANELLAR_BLAST.out.organelle_report
                                .map { it ->
                                    [[id: it[0].id, process: "MITO"], it[1]]
                                }
                                .ifEmpty { [[],[]] }

        ch_chloro           = PLASTID_ORGANELLAR_BLAST.out.organelle_report
                                .map { it ->
                                    [[id: it[0].id, process: "CHLORO"], it[1]]
                                }
                                .ifEmpty { [[],[]] }

    } else {
        ch_mito             = Channel.of( [[],[]] )
        ch_chloro           = Channel.of( [[],[]] )
    }


    //
    // SUBWORKFLOW: IDENTITY PACBIO BARCODES IN INPUT DATA
    //
    if ( (include_workflow_steps.contains('pacbio_barcodes') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("pacbio_barcodes")
    ) {
        PACBIO_BARCODE_CHECK (
            reference_tuple_from_GG,
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
    if ( (include_workflow_steps.contains('fcs-adaptor') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("fcs-adaptor")
    ) {
        RUN_FCSADAPTOR (
            reference_tuple_from_GG
        )
        ch_versions         = ch_versions.mix(RUN_FCSADAPTOR.out.versions)

        RUN_FCSADAPTOR.out.ch_euk
            .combine(
                RUN_FCSADAPTOR.out.ch_prok.map{it[1]}
            )
            .map { meta, file1, file2 ->
                tuple(
                    [id: meta.id, process: "FCS-Adaptor"],
                    file1,
                    file2
                )
            }
            .set{ ch_fcsadapt }

    } else {
        ch_fcsadapt         = Channel.of([[],[]])
    }


    //
    // SUBWORKFLOW: RUN FCS-GX TO IDENTIFY CONTAMINATION IN THE ASSEMBLY
    //
    if ( (include_workflow_steps.contains('fcs-gx') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("fcs-gx")
    ) {

        reference_tuple_from_GG
            .combine(fcs_db)
            .combine(Channel.of(params.taxid))
            .combine(Channel.of(params.ncbi_ranked_lineage_path))
            .multiMap { meta, ref, db, taxid, tax_path ->
                reference: [meta, taxid, ref]
                fcs_db_path: db
                taxid_val: taxid
                ncbi_tax_path: tax_path
            }
            .set { joint_channel }


        RUN_FCSGX (
            joint_channel.reference,
            joint_channel.fcs_db_path,
            joint_channel.ncbi_tax_path
        )
        ch_versions         = ch_versions.mix(RUN_FCSGX.out.versions)
        ch_fcsgx            = RUN_FCSGX.out.fcsgxresult
                                .map { it ->
                                    [[id: it[0].id, process: "FCSGX result"], it[1]]
                                }
                                .ifEmpty { [[],[]] }

    } else {
        ch_fcsgx         = Channel.of( [[],[]] )
    }


    //
    // SUBWORKFLOW: CALCULATE AVERAGE READ COVERAGE
    //
    if ( (include_workflow_steps.contains('coverage') || include_workflow_steps.contains('btk_busco') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("coverage")
    ) {
        RUN_READ_COVERAGE (
            reference_tuple_from_GG,
            reads,
            params.reads_type,
        )
        ch_versions         = ch_versions.mix(RUN_READ_COVERAGE.out.versions)
        ch_coverage         = RUN_READ_COVERAGE.out.tsv_ch
                                .map { it ->
                                    [[id: it[0].id, process: "Coverage"], it[1]]
                                }
                                .ifEmpty { [[],[]] }

        ch_bam              = RUN_READ_COVERAGE.out.bam_ch
                                .map { it ->
                                    [[id: it[0].id, process: "Mapped Bam"], it[1]]
                                }
                                .ifEmpty { [[],[]] }

    } else {
        ch_coverage         = Channel.of( [[],[]] )
        ch_bam              = Channel.of( [[],[]] )
    }


    //
    // SUBWORKFLOW: SCREENING FOR VECTOR SEQUENCE
    //
    if ( (include_workflow_steps.contains('vecscreen') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("vecscreen")
    ) {
        RUN_VECSCREEN (
            reference_tuple_from_GG,
            params.vecscreen_database_path
        )
        ch_versions         = ch_versions.mix(RUN_VECSCREEN.out.versions)
        ch_vecscreen        = RUN_VECSCREEN.out.vecscreen_contam
                                .map { it ->
                                    [[id: it[0].id, process: "Vecscreen"], it[1]]
                                }
                                .ifEmpty { [[],[]] }
    } else {
        ch_vecscreen        = Channel.of([[],[]])
    }


    //
    // SUBWORKFLOW: Run the kraken classifier
    //
    if ( (include_workflow_steps.contains('kraken') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("kraken")
    ) {

        RUN_NT_KRAKEN(
            reference_tuple_from_GG,
            params.nt_kraken_database_path,
            params.ncbi_ranked_lineage_path
        )
        ch_versions         = ch_versions.mix(RUN_NT_KRAKEN.out.versions)

        ch_kraken1 = RUN_NT_KRAKEN.out.classified
                        .map { it ->
                            [[id: it[0].id, process: "Kraken 1"], it[1]]
                        }
                    .ifEmpty { [[],[]] }

        ch_kraken2 = RUN_NT_KRAKEN.out.report
                        .map { it ->
                            [[id: it[0].id, process: "Kraken 2"], it[1]]
                        }
                    .ifEmpty { [[],[]] }

        ch_kraken3 = RUN_NT_KRAKEN.out.lineage
                        .map { it ->
                            [[id: it[0].id, process: "Kraken 3"], it[1]]
                        }
                    .ifEmpty { [[],[]] }
    } else {
        ch_kraken1 = Channel.of([[],[]])
        ch_kraken2 = Channel.of([[],[]])
        ch_kraken3 = Channel.of([[],[]])

    }


    if ( (include_workflow_steps.contains('create_btk_dataset') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("create_btk_dataset")
    ) {

        //
        // LOGIC: FOUND RACE CONDITION EFFECTING LONG RUNNING JOBS
        //          AND INPUT TO HERE ARE NOW MERGED AND MAPPED
        //          EMPTY CHANNELS ARE CHECKED AND DEFAULTED TO [[],[]]
        //
        reference_tuple_from_GG
            .map{ it -> tuple([
                id: it[0].id,
                taxid: it[0].taxid,
                sci_name: it[0].sci_name,
                process: "REFERENCE"], it[1])
            }
            .mix(
                ej_dot_genome.map{ it -> tuple([id: it[0].id, process: "GENOME"], it[1])},
                ch_kmers,
                ch_tiara,
                ch_nt_blast,
                // Use the BLAST top hits for BTK
                ch_btk_format,
                ch_fcsgx,
                ch_bam,
                ch_coverage,
                ch_kraken1,
                ch_kraken2,
                ch_kraken3,
                nr_full,
                un_full
            )
            .map { meta, file ->
                [meta.id, [meta: meta, file: file]]
            }
            .filter { id, data -> id != [] }
            .groupTuple()
            .map { id, data ->
                [id: id, data: data]
            }
            .set {ch_genomic_cbtk_input}


        //
        // LOGIC: LIST OF PROCESSES TO CHECK FOR
        //
        def processes = [
            'REFERENCE', 'NT-BLAST', 'TIARA', 'Kraken 2', 'GENOME', 'KMERS',
            'FCSGX result', 'NR-FULL', 'UN-FULL', 'Mapped Bam', 'Coverage',
            'Kraken 1', 'Kraken 3'
        ]


        //
        // LOGIC: Create a channel for each process
        //
        def processChannels = processes.collectEntries { process ->
            [(process): ch_genomic_cbtk_input
                .map { sample ->
                    def data = sample.data.find { it.meta.process == process }
                    data ? [sample.id, data.meta, data.file] : [sample.id, [process: process], []]
                }
            ]
        }


        //
        // LOGIC: Combine all channels using a series of combine operations
        //
        def combined_channel = processChannels['REFERENCE']
        processes.tail().each { process ->
            combined_channel = combined_channel.combine(processChannels[process], by: 0)
        }


        //
        // MODULE: CREATE A BTK COMPATIBLE DATASET FOR NEW DATA
        //
        CREATE_BTK_DATASET (
            combined_channel,
            Channel.fromPath(params.ncbi_taxonomy_path).first()
        )
        ch_versions             = ch_versions.mix(CREATE_BTK_DATASET.out.versions)
    }


    //
    // LOGIC: AUTOFILTER ASSEMBLY BY TIARA AND FCSGX RESULTS SO THE SUBWORKLOW CAN EITHER BE TRIGGERED BY THE VALUES tiara, fcs-gx, autofilter_assemlby AND EXCLUDE STEPS NOT CONTAINING autofilter_assembly
    //          OR BY include_steps CONTAINING ALL AND EXCLUDE NOT CONTAINING autofilter_assembly.
    //
    if (
        ( include_workflow_steps.contains('tiara') && include_workflow_steps.contains('fcs-gx') && include_workflow_steps.contains("autofilter_assembly") && !exclude_workflow_steps.contains("autofilter_assembly") ) ||
        ( include_workflow_steps.contains('ALL') && !exclude_workflow_steps.contains("autofilter_assembly") )
    ) {
        //
        // LOGIC: FILTER THE INPUT FOR THE AUTOFILTER STEP
        //          - We can't just combine on meta.id as some of the Channels have other data
        //              in there too so we just sanitise, and _then_ combine on 0, and
        //              _then_ add back in the taxid as we need that for this process.
        //              Thankfully taxid is a param so easy enough to add back in.
        //                  Actually, it just makes more sense to passs in as its own channel.
        //

        reference_tuple_from_GG
            .map{ it -> tuple([id: it[0].id], it[1])}
            .combine(
                ch_tiara
                    .map{ it -> tuple([id: it[0].id], it[1])},
                by: 0
            )
            .combine(
                ch_fcsgx
                    .map{ it -> tuple([id: it[0].id], it[1])},
                by: 0
            )
            .combine(
                Channel.fromPath(params.ncbi_ranked_lineage_path)
            )
            .combine(
                Channel.of(params.taxid)
            )
            .multiMap{
                meta, ref, tiara, fcs, ncbi, thetaxid ->
                    reference:  tuple([id: meta.id, taxid: thetaxid], ref)
                    tiara_file: tuple(meta, tiara)
                    fcs_file:   tuple(meta, fcs)
                    ncbi_rank:  ncbi
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
        // LOGIC: BRANCH THE CHANNEL ON WHETHER OR NOT THERE IS ABNORMAL CONTAMINATION IN THE
        //          OUTPUT FILE.
        //
        AUTOFILTER_AND_CHECK_ASSEMBLY.out.alarm_file
            .map { file -> file.text.trim() }
            .branch { it ->
                run_btk     : it.contains("YES_ABNORMAL_CONTAMINATION")
                dont_run    : true
            }
            .set { btk_bool }

            btk_bool.run_btk.view{"BTK BOOL: $it"}


        ch_versions             = ch_versions.mix(AUTOFILTER_AND_CHECK_ASSEMBLY.out.versions)
    } else {
        ch_autofilt_assem       = Channel.of([])
        ch_autofilt_indicator   = Channel.of([])
    }


    //
    // LOGIC: IF NOT IN EXCLUDE STEPS AND ( BTK_BUSCO AND AUTOFILTER IN INCLUDE STEPS _OR_ ALL STEPS ARE ACTIVE) AND MODE IS CONDITIONAL AND ABNORMAL CONTAMINATION HAS BEEN FOUND
    //              OR
    //        IF NOT IN EXCLUDE STEPS AND ( BTK_BUSCO AND AUTOFILTER IN INCLUDE STEPS _OR_ ALL STEPS ARE ACTIVE) AND MODE IS MANDATORY
    //
    if (
        (
            !exclude_workflow_steps.contains("btk_busco") &&
            ((include_workflow_steps.contains('btk_busco') && include_workflow_steps.contains("autofilter_assembly")) || include_workflow_steps.contains('ALL')) &&
            btk_busco_run_mode == "conditional" &&
            btk_bool.run_btk
        ) ||
        (
            !exclude_workflow_steps.contains("btk_busco") &&
            ((include_workflow_steps.contains('btk_busco') && include_workflow_steps.contains("autofilter_assembly") && include_workflow_steps.contains("create_btk_dataset")) || include_workflow_steps.contains('ALL')) &&
            btk_busco_run_mode == "mandatory"
        )
    ) {
        //
        // MODULE: THIS MODULE FORMATS THE INPUT DATA IN A SPECIFIC CSV FORMAT FOR
        //          USE IN THE BTK PIPELINE
        //
        GENERATE_SAMPLESHEET (
            reference_tuple_from_GG,
            params.reads_path,
            AUTOFILTER_AND_CHECK_ASSEMBLY.out.alarm_file
        )
        ch_versions         = ch_versions.mix(GENERATE_SAMPLESHEET.out.versions)

        //
        // LOGIC: STRIP THE META DATA DOWN TO id AND COMBINE ON THAT.
        //
        GENERATE_SAMPLESHEET.out.csv
            .map{ meta, csv ->
                tuple(
                    [ id: meta.id ],
                    csv
                )
            }
            .set {coverage_id}

        reference_tuple_from_GG
            .map{ it -> tuple([id:it[0].id], it[1])}
            .combine(coverage_id, by: 0)
            .multiMap { meta_1, ref, csv ->
                reference: [meta_1, ref]
                samplesheet: csv
            }
            .set { combined_input }


        //
        // PIPELINE: PREPARE THE DATA FOR USE IN THE SANGER-TOL/BLOBTOOLKIT PIPELINE
        //              WE ARE USING THE PIPELINE HERE AS A MODULE THIS REQUIRES IT
        //              TO BE USED AS A AN INTERACTIVE JOB ON WHAT EVER EXECUTOR YOU ARE USING.
        //              This will also eventually check for the above run_btk boolean from
        //              autofilter
        SANGER_TOL_BTK (
            combined_input.reference,
            combined_input.samplesheet,
            params.diamond_uniprot_database_path,
            params.nt_database_path,
            params.diamond_uniprot_database_path,
            params.ncbi_taxonomy_path,
            params.reads_path,
            params.busco_lineages_folder,
            params.busco_lineages,
            params.taxid,
        )
        ch_versions             = ch_versions.mix(SANGER_TOL_BTK.out.versions)

        //
        // LOGIC: STRIP THE META OUT OF THE REFERENCE AND CSV SO WE CAN COMBINE ON META
        //

        GENERATE_SAMPLESHEET.out.csv
            .map{ meta, file ->
                tuple([id: meta.id], file)
            }
            .set{ new_csv }


        //
        // LOGIC: COMBINE ALL THE REQUIRED CHANNELS TOGETHER INTO A MAP FOR NF-CASCADE VERSION
        //          OF SANGER_TOL_BTK TO PARSE INTO THE INPUT PARAMS
        //
        //ch_reference
        //    .combine(new_csv, by: 0)
        //    .combine(Channel.of(params.diamond_uniprot_database_path))
        //    .combine(Channel.of(params.nt_database_path))
        //    .combine(Channel.of(params.diamond_uniprot_database_path))
        //    .combine(Channel.of(params.ncbi_taxonomy_path))
        //    .combine(Channel.of(params.busco_lineages_folder))
        //    .combine(Channel.of(params.busco_lineages))
        //    .combine(Channel.of(params.taxid))
        //    .map{
        //        meta, reference, samplesheet, prot_db, nt_db, x_db, ncbi_taxdump, busco_folder, busco_lineage_vals, taxid ->
        //        [
        //            input: samplesheet,
        //            fasta: reference,
        //            blastp: prot_db,
        //            blastn: nt_db,
        //            blastx: x_db,
        //            taxdump: ncbi_taxdump,
        //            busco: busco_folder,
        //            busco_lineages: busco_lineage_vals,
        //            taxon: taxid,
        //            blastx_outext: "txt",
        //            use_work_dir_as_temp: true
        //
        //        ]
        //    }
        //    .set{ pipeline_input }

        //
        // PIPELINE: SANGER_TOL_BTK_CASCADE USES NF-CASCADE WRITTEN BY MAHESH
        //
        //SANGER_TOL_BTK_CASCADE(
        //    "sanger-tol/blobtoolkit",
        //    "-r 0.6.0 -profile sanger,singularity",
        //    [],
        //    pipeline_input,
        //    []
        //)

        //
        // MODULE: MERGE THE TWO BTK FORMATTED DATASETS INTO ONE DATASET FOR EASIER USE
        //
        CREATE_BTK_DATASET.out.btk_datasets
            .map { meta, file -> [meta.id, [meta, file]] }
            .join(
                SANGER_TOL_BTK.out.dataset
                    .map { meta, file ->
                        [meta.id, [meta, file]]
                })
            .map { id, ref, btk -> [ref[0], ref[1], btk[1]] }
            .set { merged_channel }

        MERGE_BTK_DATASETS (
            merged_channel
        )
        ch_versions             = ch_versions.mix(MERGE_BTK_DATASETS.out.versions)
        busco_merge_btk         = MERGE_BTK_DATASETS.out.busco_summary_tsv
    } else {
        busco_merge_btk         = Channel.of([[],[]])
    }


    //
    // LOGIC: EACH SUBWORKFLOW OUTPUTS EITHER AN EMPTY CHANNEL OR A FILE CHANNEL DEPENDING ON THE RUN RULES
    //          SO THE RULES FOR THIS ONLY NEED TO BE A SIMPLE "DO YOU WANT IT OR NOT"
    //
    if (
        !exclude_workflow_steps.contains("essentials") && !exclude_workflow_steps.contains("merge")
    ) {

        // Define empty channels as fallbacks for optional processes
        ch_create_btk_summary = Channel.of([[],[]])
        
        // Only use the actual output if the process was run
        if ((include_workflow_steps.contains('create_btk_dataset') || include_workflow_steps.contains('ALL')) &&
                !exclude_workflow_steps.contains("create_btk_dataset")) {
            ch_create_btk_summary = CREATE_BTK_DATASET.out.create_summary
        }

        //
        // LOGIC: FOUND RACE CONDITION EFFECTING LONG RUNNING JOBS
        //          AND INPUT TO HERE ARE NOW MERGED AND MAPPED
        //          EMPTY CHANNELS ARE CHECKED AND DEFAULTED TO [[],[]]
        //
        ej_gc_coverage
            .map{ it -> tuple([
                id: it[0].id,
                process: "GC_COV"], it[1])
            }
            .mix(
                ej_dot_genome.map{ it -> tuple([id: it[0].id, process: "GENOME"], it[1])},
                ch_create_btk_summary.map{ it -> tuple([id: it[0].id, process: "C_BTK_SUM"], it[1])},
                busco_merge_btk.map{ it -> tuple([id: it[0].id, process: "BUSCO_MERGE"], it[1])},
                ch_kmers,
                ch_tiara,
                ch_fcsgx,
                ch_coverage,
                ch_kraken3,
                ch_blast_lineage,
                nr_hits,
                un_hits
            )
            .map { meta, file ->
                [meta.id, [meta: meta, file: file]]
            }
            .filter { id, data -> id != [] }
            .groupTuple()
            .map { id, data ->
                [id: id, data: data]
            }
            .set {ascc_merged_data}

        def processes = [
            'GC_COV', 'Coverage', 'TIARA',
            'Kraken 3', 'NT-BLAST-LINEAGE', 'KMERS', 'NR-HITS', 'UN-HITS',
            'C_BTK_SUM', 'BUSCO_MERGE','FCSGX result'
        ]

        def processChannels = processes.collectEntries { process ->
            [(process): ascc_merged_data
                .map { sample ->
                    def data = sample.data.find { it.meta.process == process }
                    data ? [sample.id, data.meta, data.file] : [sample.id, [process: process], []]
                }
            ]
        }

        def ascc_combined_channels = processChannels['GC_COV']
        processes.tail().each { process ->
            ascc_combined_channels = ascc_combined_channels
                                    .combine(processChannels[process], by: 0)
        }

        //
        // SUBWORKFLOW: MERGES DATA THAT IS NOT USED IN THE CREATION OF THE BTK_DATASETS FOLDER
        //
        ASCC_MERGE_TABLES (
            ascc_combined_channels.map { it[1..-1] } // Remove the first item in tuple (mapping key)
        )
        ch_versions             = ch_versions.mix(ASCC_MERGE_TABLES.out.versions)
    }

    //
    // SUBWORKFLOW: GENERATE HTML REPORT
    //
    if ( (include_workflow_steps.contains('html_report') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("html_report")
    ) {
        // Create channels for HTML report inputs
        // Initialize all channels as empty by default
        ch_barcode_results = Channel.of([[id: "chlamydomonas_dataset_PRIMARY"],[]])
        ch_fcs_adaptor_euk = Channel.of([[id: "chlamydomonas_dataset_PRIMARY"],[]])
        ch_fcs_adaptor_prok = Channel.of([[id: "chlamydomonas_dataset_PRIMARY"],[]])
        ch_trim_ns_results = Channel.of([[id: "chlamydomonas_dataset_PRIMARY"],[]])
        ch_vecscreen_results = Channel.of([[id: "chlamydomonas_dataset_PRIMARY"],[]])
        ch_autofilter_results = Channel.of([[id: "chlamydomonas_dataset_PRIMARY"],[]])
        ch_merged_table = Channel.of([[id: "chlamydomonas_dataset_PRIMARY"],[]])
        ch_fasta_sanitation_log = Channel.of([[id: "chlamydomonas_dataset_PRIMARY"],[]])
        
        // Only access workflow outputs if the workflow was actually run
        if ((include_workflow_steps.contains('pacbio_barcodes') || include_workflow_steps.contains('ALL')) &&
                !exclude_workflow_steps.contains("pacbio_barcodes")) {
            ch_barcode_results = PACBIO_BARCODE_CHECK.out.filtered
        }
        
        if ((include_workflow_steps.contains('fcs-adaptor') || include_workflow_steps.contains('ALL')) &&
                !exclude_workflow_steps.contains("fcs-adaptor")) {
            ch_fcs_adaptor_euk = RUN_FCSADAPTOR.out.ch_euk
            ch_fcs_adaptor_prok = RUN_FCSADAPTOR.out.ch_prok
        }
        
        if (!exclude_workflow_steps.contains("essentials")) {
            ch_trim_ns_results = ESSENTIAL_JOBS.out.trailing_ns_report
        }
        
        if ((include_workflow_steps.contains('vecscreen') || include_workflow_steps.contains('ALL')) &&
                !exclude_workflow_steps.contains("vecscreen")) {
            ch_vecscreen_results = RUN_VECSCREEN.out.vecscreen_contam
        }
        
        if (((include_workflow_steps.contains('tiara') && include_workflow_steps.contains('fcs-gx') && 
              include_workflow_steps.contains("autofilter_assembly") && 
              !exclude_workflow_steps.contains("autofilter_assembly")) ||
             (include_workflow_steps.contains('ALL') && 
              !exclude_workflow_steps.contains("autofilter_assembly")))) {
            ch_autofilter_results = AUTOFILTER_AND_CHECK_ASSEMBLY.out.fcs_tiara_summary
        }
        
        if (!exclude_workflow_steps.contains("essentials") && !exclude_workflow_steps.contains("merge")) {
            ch_merged_table = ASCC_MERGE_TABLES.out.merged_table
        }
        
        // Get the FASTA sanitation log if the filter_fasta process was run
        if (!exclude_workflow_steps.contains("essentials")) {
            ch_fasta_sanitation_log = ESSENTIAL_JOBS.out.filter_fasta_sanitation_log
            ch_fasta_length_filtering_log = ESSENTIAL_JOBS.out.filter_fasta_length_filtering_log
        } else {
            ch_fasta_length_filtering_log = Channel.of([[id: "chlamydomonas_dataset_PRIMARY"],[]])
        }
        
        // Create channels for the input samplesheet and YAML parameters file
        ch_samplesheet = Channel.fromPath(params.input)
        
        // Handle the params_file parameter - check if it exists and create a channel
        if (params.containsKey('params_file') && params.params_file) {
            ch_params_file = Channel.fromPath(params.params_file)
        } else {
            ch_params_file = Channel.value([])
        }
        
        // Add debug logging
        log.info "HTML Report Generation - Input Channels in main workflow:"
        log.info "ch_barcode_results: ${ch_barcode_results.dump()}"
        log.info "ch_fcs_adaptor_euk: ${ch_fcs_adaptor_euk.dump()}"
        log.info "ch_fcs_adaptor_prok: ${ch_fcs_adaptor_prok.dump()}"
        log.info "ch_trim_ns_results: ${ch_trim_ns_results.dump()}"
        log.info "ch_vecscreen_results: ${ch_vecscreen_results.dump()}"
        log.info "ch_autofilter_results: ${ch_autofilter_results.dump()}"
        log.info "ch_merged_table: ${ch_merged_table.dump()}"
        log.info "ch_fasta_sanitation_log: ${ch_fasta_sanitation_log.dump()}"
        log.info "ch_fasta_length_filtering_log: ${ch_fasta_length_filtering_log.dump()}"
        log.info "ch_samplesheet: ${ch_samplesheet.dump()}"
        log.info "ch_params_file: ${ch_params_file.dump()}"
        
        // Get the Jinja template
        ch_jinja_template = Channel.fromPath("${baseDir}/assets/templates/ascc_report.html.jinja")

        // Create a channel with the reference file
        ch_reference_file = reference_tuple_from_GG

        GENERATE_HTML_REPORT_WORKFLOW (
            ch_barcode_results,
            ch_fcs_adaptor_euk,
            ch_fcs_adaptor_prok,
            ch_trim_ns_results,
            ch_vecscreen_results,
            ch_autofilter_results,
            ch_merged_table,
            GET_KMERS_PROFILE.out.kmers_results,
            ch_reference_file,
            ch_fasta_sanitation_log,
            ch_fasta_length_filtering_log,
            ch_jinja_template,
            ch_samplesheet,
            ch_params_file
        )
        ch_versions = ch_versions.mix(GENERATE_HTML_REPORT_WORKFLOW.out.versions)
    }

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
