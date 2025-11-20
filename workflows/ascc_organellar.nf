/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CREATE_BTK_DATASET                            } from '../modules/local/blobtoolkit/create_dataset/main'
include { AUTOFILTER_AND_CHECK_ASSEMBLY                 } from '../modules/local/autofilter/autofilter/main'

include { TIARA_TIARA                                   } from '../modules/nf-core/tiara/tiara/main'

include { ESSENTIAL_JOBS                                } from '../subworkflows/local/essential_jobs/main'
include { EXTRACT_NT_BLAST                              } from '../subworkflows/local/extract_nt_blast/main'
include { PACBIO_BARCODE_CHECK                          } from '../subworkflows/local/pacbio_barcode_check/main'
include { RUN_READ_COVERAGE                             } from '../subworkflows/local/run_read_coverage/main'
include { RUN_VECSCREEN                                 } from '../subworkflows/local/run_vecscreen/main'
include { RUN_NT_KRAKEN                                 } from '../subworkflows/local/run_nt_kraken/main'
include { RUN_FCSGX                                     } from '../subworkflows/local/run_fcsgx/main'
include { RUN_FCSADAPTOR                                } from '../subworkflows/local/run_fcsadaptor/main'
include { RUN_DIAMOND as NR_DIAMOND                     } from '../subworkflows/local/run_diamond/main'
include { RUN_DIAMOND as UP_DIAMOND                     } from '../subworkflows/local/run_diamond/main'
include { ASCC_MERGE_TABLES                             } from '../modules/local/ascc/merge_tables/main'
include { GENERATE_HTML_REPORT_WORKFLOW                 } from '../subworkflows/local/generate_html_report/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASCC_ORGANELLAR {

    take:
    ch_samplesheet          // channel: samplesheet read in from --input
    fcs_ov                  // params.fcs_override
    fcs_samplesheet         // The FCS override samplesheet for override
    fcs_db                  // [path(path)]
    reads
    scientific_name         // val(name)
    pacbio_database         // tuple [[meta.id], pacbio_database]
    ncbi_taxonomy_path
    ncbi_ranked_lineage_path
    nt_database_path
    diamond_nr_db_path
    diamond_uniprot_db_path
    taxid
    nt_kraken_db_path
    vecscreen_database_path
    reads_path
    reads_type
    ch_barcodes

    main:
    ch_versions     = channel.empty()

    //
    // LOGIC: SET run_conditional AS THE TYPICAL REQUIREMENTS TO RUN PROCESS
    //
    run_conditional = ["both", "organellar"]


    //
    // LOGIC: PRETTY NOTIFICATION OF FILES AT STAGE
    //
    ch_samplesheet
        .map { meta, sample ->
            log.info "[ASCC INFO]: ORGANELLAR WORKFLOW:\n\t-- $meta\n\t-- $sample\n"
        }

    if ( !params.run_essentials in run_conditional ) {
        log.warn("[ASCC WARN]: MAKE SURE YOU ARE AWARE YOU ARE SKIPPING ESSENTIAL JOBS, THIS INCLUDES BREAKING SCAFFOLDS OVER 1.9GB, FILTERING N\'s AND GC CONTENT REPORT (THIS WILL BREAK OTHER PROCESSES AND SHOULD ONLY BE RUN WITH `--run_essentials {both,genomic,organellar,off}`)")
    }


    //
    // SUBWORKFLOW: RUNS FILTER_FASTA, GENERATE .GENOME, CALCS GC_CONTENT AND FINDS RUNS OF N's
    //
    ESSENTIAL_JOBS(
        ch_samplesheet.filter { meta, file -> params.run_essentials in run_conditional }
    )
    ch_versions             = ch_versions.mix(ESSENTIAL_JOBS.out.versions)

    ej_reference_tuple      = ESSENTIAL_JOBS.out.reference_tuple_from_GG
                                .ifEmpty{ ch_samplesheet }
                                .map { meta, _fasta ->
                                    [[  id      : meta.id,
                                        process : "REFERENCE",
                                        sliding : params.seqkit_sliding,
                                        window  : params.seqkit_window,
                                        taxid   : params.taxid
                                    ], _fasta ]
                                }

    ej_seqkit_reference     = ESSENTIAL_JOBS.out.reference_with_seqkit
                                .ifEmpty{ ch_samplesheet }

    ej_dot_genome           = ESSENTIAL_JOBS.out.dot_genome
                                .map{ meta, _file ->
                                    [[id: meta.id, process: "GENOME"], _file]
                                }

    ej_gc_coverage          = ESSENTIAL_JOBS.out.gc_content_txt
                                .ifEmpty{ [[:],[]] }

    ej_trailing_ns          = ESSENTIAL_JOBS.out.trailing_ns_report
                                .map { meta, _file ->
                                    def new_meta = meta + [process: "TRAILING_NS"]
                                    [new_meta, _file]
                                }
                                .ifEmpty{ [[process: "TRAILING_NS"],[]] }

    ej_fasta_sanitation_log = ESSENTIAL_JOBS.out.filter_fasta_sanitation_log
                                .ifEmpty{ [[process: "REFERENCE_SANI_LOG"],[]] }

    ej_fasta_filter_log     = ESSENTIAL_JOBS.out.filter_fasta_length_filtering_log
                                .ifEmpty{ [[process: "REFERENCE_FILT_LOG"],[]] }


    // ----------------------------------------------
    //
    // SUBWORKFLOW: EXTRACT RESULTS HITS FROM TIARA
    //
    TIARA_TIARA (
        ej_reference_tuple.filter { meta, file -> params.run_tiara in run_conditional }
    )
    ch_versions         = ch_versions.mix( TIARA_TIARA.out.versions )

    ch_tiara            = TIARA_TIARA.out.classifications
                            .map { meta, file ->
                                [[id: meta.id, process: "TIARA"], file]
                            }
                            .ifEmpty{ [[process: "TIARA"],[]] }


    // ----------------------------------------------
    //
    // SUBWORKFLOW: IDENTITY PACBIO BARCODES IN INPUT DATA
    //

    ej_reference_tuple
        .filter { meta, file ->
            params.run_pacbio_barcodes in run_conditional
        }
        .combine(pacbio_database)
        .multiMap{
            ref_meta, ref_data, pdb_meta, pdb_data ->
                reference: [ref_meta, ref_data]
                pacbio_db: [pdb_meta, pdb_data]
        }
        .set { duplicated_db }

    PACBIO_BARCODE_CHECK (
        duplicated_db.reference,
        ch_barcodes,
        duplicated_db.pacbio_db
    )
    ch_versions         = ch_versions.mix(PACBIO_BARCODE_CHECK.out.versions)

    ch_barcode_check    = PACBIO_BARCODE_CHECK.out.filtered
                            .map{ meta, _file ->
                                def new_meta = meta + [process: "BARCODES"]
                                [new_meta, _file]
                            }
                            .ifEmpty{ [[process: "BARCODES"],[]] }


    // ----------------------------------------------
    //
    // SUBWORKFLOW: RUN FCS-ADAPTOR TO IDENTIDY ADAPTOR AND VECTORR CONTAMINATION
    //
    RUN_FCSADAPTOR (
        ej_reference_tuple.filter { meta, file -> params.run_fcs_adaptor in run_conditional}
    )
    ch_versions         = ch_versions.mix(RUN_FCSADAPTOR.out.versions)

    //
    // LOGIC: AT THIS POINT THE META CONTAINS JUNK THAT CAN 'CONTAMINATE' MATCHES,
    //          SO STRIP IT DOWN BEFORE USE, WE ALSO MERGE THE OUTPUT TOGETHER FOR SIMPLICITY
    //
    RUN_FCSADAPTOR.out.ch_euk
        .combine(
            RUN_FCSADAPTOR.out.ch_prok.map{it[1]}
        )
        .map { meta, file1, file2 ->
            [[ id: meta.id, process: "FCS-Adaptor"], file1, file2 ]
        }
        .ifEmpty{ [[process: "FCS-Adaptor"],[], []] }
        .set { ch_fcsadapt }


    // ----------------------------------------------
    //
    // SUBWORKFLOW: RUN FCS-GX TO IDENTIFY CONTAMINATION IN THE ASSEMBLY
    //
    if ( !params.fcs_override ) {

        joint_channel = ej_reference_tuple
            .filter { meta, file ->
                params.run_fcsgx in run_conditional
            }
            .combine(fcs_db)
            .combine(taxid)
            .combine(ncbi_ranked_lineage_path)
            .multiMap { meta, ref, db, tax_id, tax_path ->
                def new_meta = [id: meta.id, taxid: meta.taxid]
                reference: [new_meta, ref]
                fcs_db_path: db
                ncbi_tax_path: tax_path
            }

        RUN_FCSGX (
            joint_channel.reference,
            joint_channel.fcs_db_path,
            joint_channel.ncbi_tax_path
        )
        ch_versions         = ch_versions.mix(RUN_FCSGX.out.versions)

        ch_fcsgx            = RUN_FCSGX.out.fcsgxresult
                                .ifEmpty{ [[process: "FCSGX_RESULT"],[]] }

        ch_fcsgx_report     = RUN_FCSGX.out.fcsgx_report_txt
                                .ifEmpty{ [[process: "FCSGX_REPORT"],[]] }

        ch_fcsgx_taxonomy   = RUN_FCSGX.out.fcsgx_taxonomy_rpt
                                .ifEmpty{ [[process: "FCSGX_TAX_REPORT"],[]] }

    } else if ( params.fcs_override && params.run_fcsgx in ["both", "genomic"]) {

        //
        // TODO: THIS NEEDS TO BE OUPUT TO RESULTS TOO
        //
        fcs_samplesheet.map{ meta, file ->
            log.info("[ASCC INFO]: Overriding Internal FCSGX with ${file}")
            def new_meta = meta + [process: "FCSGX_RESULT"]
            [new_meta, file]

        }
        .set { ch_fcsgx }

        ch_fcsgx_report     = channel.of( [[process: "FCSGX_REPORT"],[]] )
        ch_fcsgx_taxonomy   = channel.of( [[process: "FCSGX_TAX_REPORT"],[]] )

    }


    // ----------------------------------------------
    //
    // SUBWORKFLOW: CALCULATE AVERAGE READ COVERAGE
    //

    RUN_READ_COVERAGE (
        ej_reference_tuple.filter { meta, file -> params.run_coverage in run_conditional },
        reads_path,
        reads_type.first(), //Subworkflow uses the param, not this value... as soon as it's in a channel it can't be used for a comparator.
    )
    ch_versions         = ch_versions.mix(RUN_READ_COVERAGE.out.versions)

    //
    // LOGIC: AT THIS POINT THE META CONTAINS JUNK THAT CAN 'CONTAMINATE' MATCHES,
    //          SO STRIP IT DOWN AND ADD PROCESS_NAME BEFORE USE
    //
    ch_coverage         = RUN_READ_COVERAGE.out.tsv_ch
                            .map { meta, file ->
                                [[id: meta.id, process: "Coverage"], file]
                            }
                            .ifEmpty{ [[:],[]] }

    ch_bam              = RUN_READ_COVERAGE.out.bam_ch
                            .map { meta, file ->
                                [[id: meta.id, process: "Mapped Bam"], file]
                            }
                            .ifEmpty{ [[:],[]] }


    // ----------------------------------------------
    //
    // SUBWORKFLOW: SCREENING FOR VECTOR SEQUENCE
    //
    RUN_VECSCREEN (
        ej_reference_tuple.filter { meta, file -> params.run_vecscreen in run_conditional},
        vecscreen_database_path.first()
    )
    ch_versions         = ch_versions.mix(RUN_VECSCREEN.out.versions)

    //
    // LOGIC: AT THIS POINT THE META CONTAINS JUNK THAT CAN 'CONTAMINATE' MATCHES,
    //          SO STRIP IT DOWN AND ADD PROCESS_NAME BEFORE USE
    //
    ch_vecscreen        = RUN_VECSCREEN.out.vecscreen_contam
                            .map { it ->
                                [[id: it[0].id, process: "VECSCREEN"], it[1]]
                            }
                            .ifEmpty{ [[process: "VECSCREEN"],[]] }


    // ----------------------------------------------
    //
    // SUBWORKFLOW: RUN THE KRAKEN CLASSIFIER
    //
    RUN_NT_KRAKEN(
        ej_reference_tuple.filter { meta, file -> params.run_kraken in run_conditional},
        nt_kraken_db_path.first(),
        ncbi_ranked_lineage_path.first()
    )
    ch_versions         = ch_versions.mix(RUN_NT_KRAKEN.out.versions)

    //
    // LOGIC: AT THIS POINT THE META CONTAINS JUNK THAT CAN 'CONTAMINATE' MATCHES,
    //          SO STRIP IT DOWN AND ADD PROCESS_NAME BEFORE USE
    //
    ch_kraken1 = RUN_NT_KRAKEN.out.classified
                    .map { it ->
                        [[id: it[0].id, process: "Kraken 1"], it[1]]
                    }
                    .ifEmpty{ [[:],[]] }

    ch_kraken2 = RUN_NT_KRAKEN.out.report
                    .map { it ->
                        [[id: it[0].id, process: "Kraken 2"], it[1]]
                    }
                    .ifEmpty{ [[:],[]] }

    ch_kraken3 = RUN_NT_KRAKEN.out.lineage
                    .map { it ->
                        [[id: it[0].id, process: "Kraken 3"], it[1]]
                    }
                    .ifEmpty{ [[:],[]] }


    // ----------------------------------------------
    //
    // LOGIC: WE NEED TO MAKE SURE THAT THE INPUT SEQUENCE IS OF AT LEAST LENGTH OF params.seqkit_window
    //
    valid_length_fasta = ej_seqkit_reference
        //
        // NOTE: Here we are using the un-filtered genome, any filtering may (accidently) cause an empty fasta
        //
        .map{ meta, file ->
            def total_length = 0
            file.eachLine { line ->
                if (line && !line.startsWith('>')) {
                    total_length += line.length()
                }
            }

            def meta2 = [
                    id: meta.id,
                    sliding: meta.sliding,
                    window: meta.window,
                    seq_count: total_length
                ]

            [meta2, file]
        }
        .filter { meta, file ->
                    meta.seq_count >= params.seqkit_window
        }

    valid_length_fasta
        .map{ meta, file ->
            log.info "[ASCC INFO]: Running BLAST (NT, DIAMOND, NR) on VALID ORGANELLE: \n\t-- ${meta.id}'s sequence ($meta.seq_count bases) is >= seqkit_window $params.seqkit_window\n"
        }


    // ----------------------------------------------
    //
    // LOGIC: THIS CONDITIONAL SHOULD EXECUTE THE PROCESS WHEN:
    //          INCLUDE STEPS ARE EITHER nt_blast AND all
    //              _AS WELL AS_
    //          EXCLUDE _NOT_ CONTAINING nt_blast AND THE valid_length_fasta IS NOT EMPTY
    //

    //
    // SUBWORKFLOW: EXTRACT RESULTS HITS FROM NT-BLAST
    //
    EXTRACT_NT_BLAST (
        valid_length_fasta.filter { meta, file -> params.run_nt_blast in run_conditional },
        nt_database_path.first(),
        ncbi_ranked_lineage_path.first()
    )
    ch_versions         = ch_versions.mix(EXTRACT_NT_BLAST.out.versions)
    ch_nt_blast         = EXTRACT_NT_BLAST.out.ch_blast_hits
                            .map { it ->
                                [[id: it[0].id, process: "NT-BLAST"], it[1]]
                            }
                            .ifEmpty{ [[:],[]] }

    ch_blast_lineage    = EXTRACT_NT_BLAST.out.ch_top_lineages
                            .map { it ->
                                [[id: it[0].id, process: "NT-BLAST-LINEAGE"], it[1]]
                            }
                            .ifEmpty{ [[:],[]] }

    ch_btk_format       = EXTRACT_NT_BLAST.out.ch_btk_format
                            .map { it ->
                                [[id: it[0].id, process: "NT-BLAST-BTK"], it[1]]
                            }
                            .ifEmpty{ [[:],[]] }


    // ----------------------------------------------
    //
    // SUBWORKFLOW: DIAMOND BLAST FOR INPUT ASSEMBLY
    //
    NR_DIAMOND (
        valid_length_fasta.filter { meta, file -> params.run_nr_diamond in run_conditional },
        diamond_nr_db_path.first()
    )
    ch_versions         = ch_versions.mix(NR_DIAMOND.out.versions)

    //
    // LOGIC: AT THIS POINT THE META CONTAINS JUNK THAT CAN 'CONTAMINATE' MATCHES,
    //          SO STRIP IT DOWN AND ADD PROCESS_NAME BEFORE USE
    //
    nr_full             = NR_DIAMOND.out.reformed
                            .map { it ->
                                [[id: it[0].id, process: "NR-FULL"], it[1]]
                            }
                            .ifEmpty{ [[:],[]] }

    nr_hits             = NR_DIAMOND.out.hits_file
                            .map { it ->
                                [[id: it[0].id, process: "NR-HITS"], it[1]]
                            }
                            .ifEmpty{ [[:],[]] }


    // ----------------------------------------------
    //
    // SUBWORKFLOW: DIAMOND BLAST FOR INPUT ASSEMBLY
    //
    // NOTE: Format is "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms sphylums salltitles"
    UP_DIAMOND (
        valid_length_fasta.filter { meta, file -> params.run_uniprot_diamond in run_conditional },
        diamond_uniprot_db_path.first()
    )
    ch_versions         = ch_versions.mix(UP_DIAMOND.out.versions)
    un_full             = UP_DIAMOND.out.reformed
                            .map { it ->
                                [[id: it[0].id, process: "UN-FULL"], it[1]]
                            }
                            .ifEmpty{ [[:],[]] }

    un_hits             = UP_DIAMOND.out.hits_file
                            .map { it ->
                                [[id: it[0].id, process: "UN-HITS"], it[1]]
                            }
                            .ifEmpty{ [[:],[]] }


    // ----------------------------------------------
    //
    // LOGIC: FOUND RACE CONDITION EFFECTING LONG RUNNING JOBS
    //          AND INPUT TO HERE ARE NOW MERGED AND MAPPED
    //          EMPTY CHANNELS ARE CHECKED AND DEFAULTED TO [[],[]]
    //
    ch_organellar_cbtk_input = ej_reference_tuple
        .filter { meta, file ->
            params.run_create_btk_dataset in run_conditional
        }
        .map{ meta, _file ->
            [[
                id: meta.id,
                taxid: meta.taxid,
                sci_name: meta.sci_name,
                process: "REFERENCE"
            ], _file]
        }

        .mix(
            ej_reference_tuple.map{ it -> [[id: it[0].id, process: "GENOME"], it[1]] },
            ch_tiara,
            ch_nt_blast,
            ch_btk_format,
            // ch_fcs
            // ch_kmers were removed
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
        .filter { id, data -> id != [] && id != null }
        .groupTuple()
        .map { id, data ->
            [id: id, data: data]
        }


    //
    // LOGIC: LIST OF PROCESSES TO CHECK FOR
    //
    def processes_cbd = [
        'REFERENCE', 'NT-BLAST', 'TIARA', 'Kraken 2', 'GENOME', 'KMERS',
        'FCSGX_RESULT', 'NR-FULL', 'UN-FULL', 'Mapped Bam', 'Coverage',
        'Kraken 1', 'Kraken 3'
    ]

    //
    // LOGIC: Create a channel for each process
    //
    def processChannels_cbd = processes_cbd.collectEntries { process ->
        [(process): ch_organellar_cbtk_input
            .map { sample ->
                def data = sample.data.find { it.meta.process == process }
                data ? [sample.id, data.meta, data.file] : [sample.id, [process: process], []]
            }
        ]
    }


    //
    // LOGIC: Combine all channels using a series of combine operations
    //
    def combined_channel = processChannels_cbd['REFERENCE']

    processes_cbd.tail().each { process ->
        combined_channel = combined_channel.combine(processChannels_cbd[process], by: 0)
    }


    //
    // MODULE: CREATE A BTK COMPATIBLE DATASET FOR NEW DATA
    //
    CREATE_BTK_DATASET (
        combined_channel,
        ncbi_taxonomy_path.first(),
        scientific_name
    )
    ch_versions             = ch_versions.mix(CREATE_BTK_DATASET.out.versions)

    ch_create_summary       = CREATE_BTK_DATASET.out.create_summary
                                .map{ meta, file ->
                                    [[id: meta.id, process: "C_BTK_SUM"], file]
                                }
    ch_create_btk_dataset   = CREATE_BTK_DATASET.out.btk_datasets
                                .map{ meta, _file ->
                                    def new_meta = meta + [process: "BTK_DATASET"]
                                    [new_meta, _file]
                                }


    //
    // LOGIC: AUTOFILTER ASSEMBLY BY TIARA AND FCSGX RESULTS SO THE SUBWORKLOW CAN EITHER BE TRIGGERED BY THE VALUES tiara, fcs-gx, autofilter_assemlby AND EXCLUDE STEPS NOT CONTAINING autofilter_assembly
    //          OR BY include_steps CONTAINING ALL AND EXCLUDE NOT CONTAINING autofilter_assembly.
    //
    if (
        ( params.run_tiara == "both" || params.run_tiara == "organellar" ) &&
        ( params.run_fcsgx == "both" || params.run_fcsgx == "organellar" ) &&
        ( params.run_autofilter_assembly == "both" || params.run_autofilter_assembly == "organellar" )
    ) {
        //
        // LOGIC: FILTER THE INPUT FOR THE AUTOFILTER STEP
        //          - We can't just combine on meta.id as some of the Channels have other data
        //              in there too so we just sanitise, and _then_ combine on 0, and
        //              _then_ add back in the taxid as we need that for this process.
        //              Thankfully taxid is a param so easy enough to add back in.
        //                  Actually, it just makes more sense to passs in as its own channel.
        //

        autofilter_input_formatted = ej_reference_tuple
            .map{ it -> [[id: it[0].id], it[1]] }
            .combine(
                ch_tiara
                    .map{ it -> [[id: it[0].id], it[1]] },
                by: 0
            )
            .combine(
                ch_fcsgx
                    .map{ it -> [[id: it[0].id], it[1]] },
                by: 0
            )
            .combine(
                ncbi_ranked_lineage_path
            )
            .combine(
                taxid
            )
            .multiMap{
                meta, ref, tiara, fcs, ncbi, thetaxid ->
                    def new_meta = [id: meta.id, taxid: thetaxid]
                    reference:  [new_meta, ref]
                    tiara_file: [new_meta, tiara]
                    fcs_file:   [new_meta, fcs]
                    ncbi_rank:  ncbi
            }


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

        ch_autofilt_alarm_file  = AUTOFILTER_AND_CHECK_ASSEMBLY.out.alarm_file
            .map{ meta, file ->
                [[id: meta.id], file]
            }

        ch_autofilt_fcs_tiara   = AUTOFILTER_AND_CHECK_ASSEMBLY.out.fcs_tiara_summary
                                    .map {meta, _file ->
                                        def new_meta = meta + [process: "AUTOFILTER"]
                                        [new_meta, _file]
                                    }
        ch_autofilt_removed_seqs= AUTOFILTER_AND_CHECK_ASSEMBLY.out.removed_seqs
        ch_autofilt_raw_report  = AUTOFILTER_AND_CHECK_ASSEMBLY.out.raw_report

        ch_versions             = ch_versions.mix(AUTOFILTER_AND_CHECK_ASSEMBLY.out.versions)
    } else {
        ch_autofilt_alarm_file  = channel.of( [[:],[]] )
        ch_autofilt_removed_seqs= channel.of( [[:],[]] )
        ch_autofilt_assem       = channel.of( [[:],[]] )
        ch_autofilt_indicator   = channel.of( [[:],[]] )
        ch_autofilt_fcs_tiara   = channel.of( [[process: "AUTOFILTER"],[]] )
        ch_autofilt_raw_report  = channel.of( [[:],[]] )
    }

    //
    // LOGIC: EACH SUBWORKFLOW OUTPUTS EITHER AN EMPTY CHANNEL OR A FILE CHANNEL DEPENDING ON THE RUN RULES
    //          SO THE RULES FOR THIS ONLY NEED TO BE A SIMPLE "DO YOU WANT IT OR NOT"
    //
    if (
        ( params.run_essentials == "both" || params.run_essentials == "organellar" ) &&
        ( params.run_merge_datasets == "both" || params.run_merge_datasets == "organellar" )
    ) {

        //
        // LOGIC: FOUND RACE CONDITION EFFECTING LONG RUNNING JOBS
        //          AND INPUT TO HERE ARE NOW MERGED AND MAPPED
        //          EMPTY CHANNELS ARE CHECKED AND DEFAULTED TO [[],[]]
        //
        busco_merge_btk = channel.of( [[],[]] )
        ch_kmers        = channel.of( [[],[]] )

        ascc_merged_data = ej_gc_coverage
            .map{ meta, file ->
                [[id: meta.id, process: "GC_COV"], file]
            }
            .mix(
                ej_dot_genome,
                ch_create_summary,
                busco_merge_btk,
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

        def processes_mbd = [
            'GC_COV', 'Coverage', 'TIARA',
            'Kraken 3', 'NT-BLAST-LINEAGE', 'KMERS', 'NR-HITS', 'UN-HITS',
            'C_BTK_SUM', 'BUSCO_MERGE','FCSGX_RESULT'
        ]

        def processChannels_mbd = processes_mbd.collectEntries { process ->
            [(process): ascc_merged_data
                .map { sample ->
                    def data = sample.data.find { it.meta.process == process }
                    data ? [sample.id, data.meta, data.file] : [sample.id, [process: process], []]
                }
            ]
        }

        def ascc_combined_channels = processChannels_mbd['GC_COV']

        processes_mbd.tail().each { process ->
            ascc_combined_channels = ascc_combined_channels
                                    .combine(processChannels_mbd[process], by: 0)
        }

        ASCC_MERGE_TABLES (
            ascc_combined_channels.map { it[1..-1] }
        )
        ch_versions               = ch_versions.mix(ASCC_MERGE_TABLES.out.versions)
        org_merged_table          = ASCC_MERGE_TABLES.out.merged_table
                                        .map { meta, _file ->
                                            def new_meta = meta + [process: "MERGED_TABLE"]
                                            [new_meta, _file]
                                        }

        org_merged_phylum_count   = ASCC_MERGE_TABLES.out.phylum_counts
                                        .map { meta, _file ->
                                            [[id: meta.id, process: "MERGED_PHYLUM_COUNTS"], _file]
                                        }

    } else {
        org_merged_table          = channel.from( [[process: "MERGED_TABLE"],[]] )
        merged_extended_table     = channel.empty()
        org_merged_phylum_count   = channel.from( [[process: "MERGED_PHYLUM_COUNTS"],[]] )
    }

    // ----------------------------------------------
    //
    // SUBWORKFLOW: GENERATE HTML REPORT (minimal wiring, opt-in)
    //  Gate with params.run_html_report to avoid altering default behavior.
    //  Use existing channels; substitute placeholders where features are disabled or not wired.
    //

    // Placeholders for channels not generated in organellar subworkflow
    ch_kmers_results = channel.of( [[],[]] )

    // Samplesheet/params file
    ch_samplesheet_path = channel.fromPath(params.input)
    ch_params_file      = params.params_file ? channel.fromPath(params.params_file) : channel.value([])

    // Templates and CSS
    ch_jinja_templates = channel.fromPath("${baseDir}/assets/templates/*.jinja").collect()
    ch_css_files       = channel.fromPath("${baseDir}/assets/css/*.css").collect()

    // GENERATE_HTML_REPORT_WORKFLOW (
    //     ch_barcode_check,
    //     ch_fcsadapt,
    //     ej_trailing_ns,
    //     ch_vecscreen,
    //     ch_autofilt_fcs_tiara,
    //     org_merged_table,
    //     org_merged_phylum_count,
    //     ch_kmers_results,
    //     ej_reference_tuple.filter {meta, file -> params.run_html_report in ["both", "organellar"]},
    //     ej_fasta_sanitation_log,
    //     ej_fasta_filter_log,
    //     ch_jinja_templates,
    //     ch_samplesheet_path,
    //     ch_params_file,
    //     ch_fcsgx_report,
    //     ch_fcsgx_taxonomy,
    //     ch_create_btk_dataset,
    //     ch_css_files
    // )
    // ch_versions = ch_versions.mix(GENERATE_HTML_REPORT_WORKFLOW.out.versions)


    emit:

    essential_reference         = ej_reference_tuple
    essential_genome_file       = ej_dot_genome
    essential_gc_cov            = ej_gc_coverage

    blast_output                = ch_nt_blast
    blast_lineage               = ch_blast_lineage
    blast_btk_formatted         = ch_btk_format

    diamond_nr_blast_full       = nr_full
    diamond_nr_blast_hits       = nr_hits

    diamond_un_blast_full       = un_full
    diamond_un_blast_hits       = un_hits

    read_coverage_output        = ch_coverage
    read_coverage_bam           = ch_bam

    fcsadaptor_prok_euk         = ch_fcsadapt
    fcsgx_output                = ch_fcsgx

    autofilter_deconned_assm    = ch_autofilt_assem
    autofilter_fcs_tiar_smry    = ch_autofilt_fcs_tiara
    autofilter_removed_seqs     = ch_autofilt_removed_seqs
    autofilter_alarm_file       = ch_autofilt_alarm_file
    autofilter_indicator_file   = ch_autofilt_indicator
    autofilter_raw_report       = ch_autofilt_raw_report

    create_btk_ds_dataset       = ch_create_btk_dataset
    create_btk_ds_create_smry   = ch_create_summary

    kraken2_classified          = ch_kraken1
    kraken2_report              = ch_kraken2
    kraken2_lineage             = ch_kraken3

    vecscreen_contam            = ch_vecscreen

    tiara_output                = ch_tiara

    versions                    = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
