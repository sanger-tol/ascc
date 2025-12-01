/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CREATE_BTK_DATASET                            } from '../modules/local/blobtoolkit/create_dataset/main'
include { MERGE_BTK_DATASETS                            } from '../modules/local/blobtoolkit/merge_dataset/main'
include { ASCC_MERGE_TABLES                             } from '../modules/local/ascc/merge_tables/main'
include { AUTOFILTER_AND_CHECK_ASSEMBLY                 } from '../modules/local/autofilter/autofilter/main'
include { SANGER_TOL_BTK                                } from '../modules/local/sanger-tol/btk/main'
include { BLOBTOOLKIT_GENERATECSV                       } from '../modules/sanger-tol/blobtoolkit/generatecsv/main'

include { TIARA_TIARA                                   } from '../modules/nf-core/tiara/tiara/main'

include { ESSENTIAL_JOBS                                } from '../subworkflows/local/essential_jobs/main'
include { GET_KMERS_PROFILE                             } from '../subworkflows/local/get_kmers_profile/main'
include { EXTRACT_NT_BLAST                              } from '../subworkflows/local/extract_nt_blast/main'
include { ORGANELLAR_BLAST as PLASTID_ORGANELLAR_BLAST  } from '../subworkflows/local/organellar_blast/main'
include { ORGANELLAR_BLAST as MITO_ORGANELLAR_BLAST     } from '../subworkflows/local/organellar_blast/main'
include { PACBIO_BARCODE_CHECK                          } from '../subworkflows/local/pacbio_barcode_check/main'
include { RUN_READ_COVERAGE                             } from '../subworkflows/local/run_read_coverage/main'
include { RUN_VECSCREEN                                 } from '../subworkflows/local/run_vecscreen/main'
include { RUN_NT_KRAKEN                                 } from '../subworkflows/local/run_nt_kraken/main'
include { RUN_FCSGX                                     } from '../subworkflows/local/run_fcsgx/main'
include { RUN_FCSADAPTOR                                } from '../subworkflows/local/run_fcsadaptor/main'
include { RUN_DIAMOND as NR_DIAMOND                     } from '../subworkflows/local/run_diamond/main'
include { RUN_DIAMOND as UP_DIAMOND                     } from '../subworkflows/local/run_diamond/main'
include { GENERATE_HTML_REPORT_WORKFLOW                 } from '../subworkflows/local/generate_html_report/main'

include { getEmptyPlaceholder } from "${projectDir}/lib/ascc_utils.groovy"



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASCC_GENOMIC {

    take:
    ch_samplesheet          // channel: samplesheet read in from --input
    organellar_genomes      // channel: tuple(meta, reference)
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
    reads_layout
    reads_type
    btk_lineages
    btk_lineages_path
    ch_barcodes

    main:
    ch_versions = channel.empty()

    //
    // LOGIC: CREATE run_conditional LIST
    //
    run_conditionals = ["both", "genomic"]


    //
    // LOGIC: PRETTY NOTIFICATION OF FILES AT STAGE
    //
    ch_samplesheet
        .map { meta, sample ->
            log.info "[ASCC INFO]: GENOMIC WORKFLOW:\n\t-- $meta\n\t-- $sample\n"
        }


    //-------------------------------------------------------------------------
    //
    // SUBWORKFLOW: RUNS FILTER_FASTA, GENERATE .GENOME, CALCS GC_CONTENT AND FINDS RUNS OF N's
    //                  THIS SHOULD NOT RUN ONLY WHEN SPECIFICALLY REQUESTED
    //
    ESSENTIAL_JOBS(
        ch_samplesheet
    )
    ch_versions             = ch_versions.mix(ESSENTIAL_JOBS.out.versions)
    ej_reference_tuple      = ESSENTIAL_JOBS.out.reference_tuple_from_GG
    ej_dot_genome           = ESSENTIAL_JOBS.out.dot_genome
    ej_gc_coverage          = ESSENTIAL_JOBS.out.gc_content_txt
    ej_trailing_ns          = ESSENTIAL_JOBS.out.trailing_ns_report
    ej_fasta_sanitation_log = ESSENTIAL_JOBS.out.filter_fasta_sanitation_log
    ej_fasta_filter_log     = ESSENTIAL_JOBS.out.filter_fasta_length_filtering_log


    //-------------------------------------------------------------------------
    //
    // LOGIC: CONVERT THE CHANNEL I AN EPOCH COUNT FOR THE GET_KMER_PROFILE
    //
    ej_reference_tuple
        .map { meta, file ->
            file.countFasta() * 3
        }
        .set {autoencoder_epochs_count}

    //
    // SUBWORKFLOW: COUNT KMERS, THEN REDUCE DIMENSIONS USING SELECTED METHODS
    //
    GET_KMERS_PROFILE (
        ej_reference_tuple.filter{ meta, file -> params.run_kmers in run_conditionals },
        params.kmer_length,
        params.dimensionality_reduction_methods,
        autoencoder_epochs_count
    )
    ch_versions         = ch_versions.mix(GET_KMERS_PROFILE.out.versions)

    //
    // LOGIC: AT THIS POINT THE META CONTAINS JUNK THAT CAN 'CONTAMINATE' MATCHES,
    //          SO STRIP IT DOWN AND ADD PROCESS_NAME BEFORE USE
    //
    ch_kmers            = GET_KMERS_PROFILE.out.combined_csv.ifEmpty { [[:],[]] }
    ch_kmers_results    = GET_KMERS_PROFILE.out.kmers_results.ifEmpty { [[:],[]] }


    //-------------------------------------------------------------------------
    //
    // SUBWORKFLOW: EXTRACT RESULTS HITS FROM TIARA
    //
    TIARA_TIARA (
        ej_reference_tuple.filter{ meta, file -> params.run_tiara in run_conditionals }
    )
    ch_versions         = ch_versions.mix( TIARA_TIARA.out.versions )
    ch_tiara            = TIARA_TIARA.out.classifications
                            .map { meta, file -> [[id: meta.id ], file] }
                            .ifEmpty { [[:],[]] }


    //-------------------------------------------------------------------------
    //
    // SUBWORKFLOW: EXTRACT RESULTS HITS FROM NT-BLAST
    //

    EXTRACT_NT_BLAST (
        ej_reference_tuple.filter{ meta, file ->
            params.run_nt_blast in run_conditionals
        },
        nt_database_path.first(),
        ncbi_ranked_lineage_path.first()
    )
    ch_versions         = ch_versions.mix(EXTRACT_NT_BLAST.out.versions)
    ch_nt_blast         = EXTRACT_NT_BLAST.out.ch_blast_hits.ifEmpty { [[:],[]] }
    ch_blast_lineage    = EXTRACT_NT_BLAST.out.ch_top_lineages.ifEmpty { [[:],[]] }
    ch_btk_format       = EXTRACT_NT_BLAST.out.ch_btk_format.ifEmpty { [[:],[]] }


    //-------------------------------------------------------------------------
    //
    // SUBWORKFLOW: DIAMOND BLAST FOR INPUT ASSEMBLY
    //

    NR_DIAMOND (
        ej_reference_tuple.filter{ meta, file ->
            params.run_nr_diamond in run_conditionals
        },
        diamond_nr_db_path.first()
    )
    ch_versions = ch_versions.mix(NR_DIAMOND.out.versions)
    nr_full     = NR_DIAMOND.out.reformed
                    .map { meta, file -> [[id: meta.id ], file] }
                    .ifEmpty { [[:],[]] }

    nr_hits     = NR_DIAMOND.out.hits_file
                    .map { meta, file -> [[id: meta.id ], file] }
                    .ifEmpty { [[:],[]] }


    //-------------------------------------------------------------------------
    //
    // SUBWORKFLOW: DIAMOND BLAST FOR INPUT ASSEMBLY
    //
    // NOTE: HEADER FORMAT WILL BE -
    //  qseqid sseqid pident length mismatch gapopen qstart qend sstart send
    //  evalue bitscore staxids sscinames sskingdoms sphylums salltitles
    UP_DIAMOND (
        ej_reference_tuple.filter{ meta, file ->
            params.run_uniprot_diamond in run_conditionals
        },
        diamond_uniprot_db_path.first()
    )
    ch_versions = ch_versions.mix(UP_DIAMOND.out.versions)

    un_full     = UP_DIAMOND.out.reformed
                    .map { meta, file -> [[id: meta.id], file ] }
                    .ifEmpty { [[:],[]] }

    un_hits     = UP_DIAMOND.out.hits_file
                    .map { meta, file -> [[id: meta.id ], file ] }
                    .ifEmpty { [[:],[]] }


    //-------------------------------------------------------------------------
    //
    // LOGIC: CHECK WHETHER THERE IS A MITO AND BRANCH
    //
    organellar_check = organellar_genomes
        .filter{ meta, file ->
            params.run_organellar_blast in run_conditionals
        }
        .branch { meta, assembly ->
            mito:       meta.assembly_type == "MITO"
            plastid:    meta.assembly_type == "PLASTID"
            invalid:    true    // if value but not of the above conditions
        }


    //
    // SUBWORKFLOW: BLASTING FOR MITO ASSEMBLIES IN GENOME
    //
    MITO_ORGANELLAR_BLAST (
        ej_reference_tuple,
        organellar_check.mito
    )
    ch_versions = ch_versions.mix(MITO_ORGANELLAR_BLAST.out.versions)


    //
    // SUBWORKFLOW: BLASTING FOR PLASTID ASSEMBLIES IN GENOME
    //
    PLASTID_ORGANELLAR_BLAST (
        ej_reference_tuple,
        organellar_check.plastid
    )
    ch_versions = ch_versions.mix(PLASTID_ORGANELLAR_BLAST.out.versions)


    //
    // LOGIC: AT THIS POINT THE META CONTAINS JUNK THAT CAN 'CONTAMINATE' MATCHES,
    //          SO STRIP IT DOWN AND ADD PROCESS_NAME BEFORE USE
    //
    ch_mito     = MITO_ORGANELLAR_BLAST.out.organelle_report
                    .map { meta, file -> [[id: meta.id ], file] }
                    .ifEmpty { [[:],[]] }

    ch_chloro   = PLASTID_ORGANELLAR_BLAST.out.organelle_report
                    .map { meta, file -> [[id: meta.id ], file] }
                    .ifEmpty { [[:],[]] }



    //-------------------------------------------------------------------------
    //
    // SUBWORKFLOW: IDENTITY PACBIO BARCODES IN INPUT DATA
    //
    ej_reference_tuple
        .combine(pacbio_database)
        .multiMap{
            ref_meta, ref_data, pdb_meta, pdb_data ->
                reference: [ref_meta, ref_data]
                pacbio_db: [pdb_meta, pdb_data]
        }
        .set { duplicated_db }

    PACBIO_BARCODE_CHECK (
        duplicated_db.reference.filter{ meta, file ->
            params.run_pacbio_barcodes in run_conditionals
        },
        ch_barcodes,
        duplicated_db.pacbio_db
    )
    ch_versions         = ch_versions.mix(PACBIO_BARCODE_CHECK.out.versions)
    ch_barcode_check    = PACBIO_BARCODE_CHECK.out.filtered.ifEmpty{ [[:],[]] }


    //-------------------------------------------------------------------------
    //
    // SUBWORKFLOW: RUN FCS-ADAPTOR TO IDENTIDY ADAPTOR AND VECTORR CONTAMINATION
    //
    RUN_FCSADAPTOR (
        ej_reference_tuple.filter{ meta, file ->
            params.run_fcs_adaptor in run_conditionals
        }
    )
    ch_versions = ch_versions.mix(RUN_FCSADAPTOR.out.versions)
    ch_fcsadapt = RUN_FCSADAPTOR.out.ch_joint_report.ifEmpty{ [[:],[]] }


    //-------------------------------------------------------------------------
    //
    // SUBWORKFLOW: RUN FCS-GX TO IDENTIFY CONTAMINATION IN THE ASSEMBLY
    //
    if ( params.run_fcsgx in run_conditionals && !params.fcs_override ) {

        joint_channel = ej_reference_tuple
            .combine(fcs_db)
            .combine(taxid)
            .combine(ncbi_ranked_lineage_path)
            .multiMap { meta, ref, db, tax_id, tax_path ->
                new_meta = [id: meta.id, taxid: meta.taxid]
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
        ch_fcsgx_report     = RUN_FCSGX.out.fcsgx_report_txt
        ch_fcsgx_taxonomy   = RUN_FCSGX.out.fcsgx_taxonomy_rpt

    } else if ( params.fcs_override ) {

        fcs_samplesheet.map{ meta, file ->
            log.info("[ASCC INFO]: Overriding Internal FCSGX with ${file}")
            [[id: meta.id], file]

        }
        .set { ch_fcsgx }

        ch_fcsgx_report     = channel.of( [[:],[]] )
        ch_fcsgx_taxonomy   = channel.of( [[:],[]] )

    } else {
        ch_fcsgx            = channel.of( [[:],[]] )
        ch_fcsgx_report     = channel.of( [[:],[]] )
        ch_fcsgx_taxonomy   = channel.of( [[:],[]] )
    }


    //-------------------------------------------------------------------------
    //
    // SUBWORKFLOW: CALCULATE AVERAGE READ COVERAGE
    //
    RUN_READ_COVERAGE (
        ej_reference_tuple.filter{ meta, file ->
            params.run_coverage in run_conditionals
        },
        reads_path,
        reads_type.first(), //Subworkflow uses the param, not this value... as soon as it's in a channel it can't be used for a comparator.
    )
    ch_versions         = ch_versions.mix(RUN_READ_COVERAGE.out.versions)
    ch_coverage         = RUN_READ_COVERAGE.out.tsv_ch.ifEmpty{ [[:], []] }
    ch_bam              = RUN_READ_COVERAGE.out.bam_ch.ifEmpty{ [[:], []] }


    //-------------------------------------------------------------------------
    //
    // SUBWORKFLOW: SCREENING FOR VECTOR SEQUENCE
    //
    RUN_VECSCREEN (
        ej_reference_tuple.filter{ meta, file ->
            params.run_vecscreen in run_conditionals
        },
        vecscreen_database_path.first()
    )
    ch_versions         = ch_versions.mix(RUN_VECSCREEN.out.versions)
    ch_vecscreen        = RUN_VECSCREEN.out.vecscreen_contam.ifEmpty{ [[:],[]] }


    //-------------------------------------------------------------------------
    //
    // SUBWORKFLOW: RUN THE KRAKEN CLASSIFIER
    //
    RUN_NT_KRAKEN(
        ej_reference_tuple.filter{ meta, file ->
            params.run_kraken in run_conditionals
        },
        nt_kraken_db_path.first(),
        ncbi_ranked_lineage_path.first()
    )
    ch_versions         = ch_versions.mix(RUN_NT_KRAKEN.out.versions)
    ch_kraken1          = RUN_NT_KRAKEN.out.classified.ifEmpty{ [[:], []] }
    ch_kraken2          = RUN_NT_KRAKEN.out.report.ifEmpty{ [[:], []] }
    ch_kraken3          = RUN_NT_KRAKEN.out.lineage.ifEmpty{ [[:], []] }


    //-------------------------------------------------------------------------
    if ( params.run_create_btk_dataset in run_conditionals ) {

        //
        // LOGIC: FOUND RACE CONDITION EFFECTING LONG RUNNING JOBS
        //          AND INPUT TO HERE ARE NOW MERGED AND MAPPED
        //          EMPTY CHANNELS ARE CHECKED AND DEFAULTED TO [[:],[]]
        //
        //
        ej_reference_tuple
            .map{meta, file -> [[id: meta.id], file]}
            .join(ch_nt_blast,  remainder: true)
            .join(ch_tiara,     remainder: true)
            .join(ej_dot_genome,remainder: true)
            .join(ch_fcsgx,     remainder: true)
            .join(ch_bam,       remainder: true)
            .join(ch_coverage,  remainder: true)
            .join(ch_kmers,     remainder: true)
            .join(ch_kraken1,   remainder: true)
            .join(ch_kraken2,   remainder: true)
            .join(ch_kraken3,   remainder: true)
            .join(nr_full,      remainder: true)
            .join(un_full,      remainder: true)
            .filter { items ->
                def meta = items[0]
                meta != null &&
                meta != [] &&
                !(meta instanceof Map && (meta.id == null || meta.isEmpty()))
            }
            .map { items ->
                // Replace null values with placeholder file
                items.withIndex().collect { item, index ->
                    if (item == null) {
                        getEmptyPlaceholder(index)
                    } else if (item instanceof List && item.isEmpty()) {
                        getEmptyPlaceholder(index)
                    } else {
                        item
                    }
                }
            }
            .set{ create_input_channel}


        //
        // MODULE: CREATE A BTK COMPATIBLE DATASET FOR NEW DATA
        //
        CREATE_BTK_DATASET (
            create_input_channel,
            params.taxid,
            ncbi_taxonomy_path.first(),
            scientific_name

        )
        ch_versions             = ch_versions.mix(CREATE_BTK_DATASET.out.versions)

        ch_create_summary       = CREATE_BTK_DATASET.out.create_summary
                                    .map{ meta, _file -> [[ id: meta.id ], _file] }

        ch_create_btk_dataset   = CREATE_BTK_DATASET.out.btk_datasets
                                    .map{ meta, _file -> [[ id: meta.id ], _file] }
    } else {
        ch_create_summary       = channel.of( [[:],[]] )
        ch_create_btk_dataset   = channel.of( [[:],[]] )
    }


    //-------------------------------------------------------------------------
    //
    // LOGIC: AUTOFILTER ASSEMBLY BY TIARA AND FCSGX RESULTS SO THE SUBWORKLOW CAN EITHER BE TRIGGERED BY THE VALUES tiara, fcs-gx, autofilter_assemlby AND EXCLUDE STEPS NOT CONTAINING autofilter_assembly
    //          OR BY include_steps CONTAINING ALL AND EXCLUDE NOT CONTAINING autofilter_assembly.
    //
    if (
        ( params.run_tiara in run_conditionals ) &&
        ( params.run_fcsgx in run_conditionals ) &&
        ( params.run_autofilter_assembly in run_conditionals )
    ) {
        //
        // LOGIC: FILTER THE INPUT FOR THE AUTOFILTER STEP
        //          - We can't just combine on meta.id as some of the channel. have other data
        //              in there too so we just sanitise, and _then_ combine on 0, and
        //              _then_ add back in the taxid as we need that for this process.
        //              Thankfully taxid is a param so easy enough to add back in.
        //                  Actually, it just makes more sense to passs in as its own channel.
        //

        ej_reference_tuple
            .map{ meta, file -> [[id: meta.id], file] }
            .combine(
                ch_tiara.map{ meta, file -> [[id: meta.id], file] },
                by: 0
            )
            .combine(
                ch_fcsgx.map{ meta, file -> [[id: meta.id], file] },
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
            .set { autofilter_input_formatted }

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
        //          CHANGE OUTPUT NAME TO BE REFERENCE NAME AND THEN ALARM FILE
        btk_bool = AUTOFILTER_AND_CHECK_ASSEMBLY.out.alarm_file
            .map { meta, file -> [meta, file.text.trim()] }
            .branch { meta, data ->
                log.info("[ASCC INFO]: Run for ${meta.id} has:\n${data}\n")

                run_btk     : data.contains("YES_ABNORMAL_CONTAMINATION")
                dont_run    : true // only other lines to be produced are "NO_ABNORMAL_CONTAMINATION"
            }
        ch_versions             = ch_versions.mix(AUTOFILTER_AND_CHECK_ASSEMBLY.out.versions)
        btk_bool_run_btk        = btk_bool.run_btk
        ch_autofilt_removed_seqs= AUTOFILTER_AND_CHECK_ASSEMBLY.out.removed_seqs
        ch_autofilt_raw_report  = AUTOFILTER_AND_CHECK_ASSEMBLY.out.raw_report

        ch_autofilt_alarm_file  = AUTOFILTER_AND_CHECK_ASSEMBLY.out.alarm_file
                                    .map{ meta, file -> [[id: meta.id], file ] }

        ch_autofilt_fcs_tiara   = AUTOFILTER_AND_CHECK_ASSEMBLY.out.fcs_tiara_summary
                                    .map{ meta, _file -> [[id: meta.id], _file] }

    } else {
        btk_bool_run_btk        = channel.of([[id: "NA"], "false"])
        ch_autofilt_alarm_file  = channel.of( [[:],[]] )
        ch_autofilt_removed_seqs= channel.of( [[:],[]] )
        ch_autofilt_assem       = channel.of( [[:],[]] )
        ch_autofilt_indicator   = channel.of( [[:],[]] )
        ch_autofilt_fcs_tiara   = channel.of( [[:],[]] )
        ch_autofilt_raw_report  = channel.of( [[:],[]] )
    }


    //-------------------------------------------------------------------------
    //
    // LOGIC: DETERMINE WHETHER BLOBTOOLKIT SHOULD RUN BASED ON CONDITIONALS
    //         - ALWAYS RUN IF params.btk_busco_run_mode == "mandatory" AND BTK

    run_btk_conditional = ej_reference_tuple
        .map { meta, file ->
                [[id: meta.id, taxid: meta.taxid], file]
            }
        // below is combined into the tuple to enforce the block to only run when channel is present.
        .combine ( btk_bool_run_btk
                        .map{ meta, data ->
                            def joined_content = data
                                    .replaceAll(/\s*\|\s*/, "-")     // Replace " | " with "-"
                                    .replaceAll(/\s+/, "-")          // Replace remaining spaces with "-"
                                    .replaceAll(/_+/, "_")           // Keep underscores as they are
                                    .replaceAll(/-+/, "-")           // Clean up multiple dashes
                            [[id: meta.id, taxid: meta.taxid], joined_content]
                        },
                by: [0]
            )
        .branch { meta, assembly, data ->
            def btk_requested           = params.run_btk_busco == "both" || params.run_btk_busco == "genomic"
            def autofilter_requested    = params.run_autofilter_assembly == "both" || params.run_autofilter_assembly == "genomic"

            def ignore_autofilter       = params.btk_busco_run_mode == "mandatory" && btk_requested
            def not_mandatory_btk       = params.btk_busco_run_mode == "conditional" && autofilter_requested && btk_requested && data.contains("YES_ABNORMAL_CONTAMINATION")

            run_btk: (ignore_autofilter || not_mandatory_btk)
            skip_btk: true
        }

    run_btk_conditional.skip_btk
        .map { meta, file, data ->
            log.info "[ASCC INFO]: CONTAMINATION THRESHOLD NOT MET"
            log.info "\t- SKIPPING BLOBTOOLKIT FOR: $meta.id"
            log.info "\t- You can verify here: $file"
        }

    if (params.run_autofilter_assembly == "off" && params.run_btk_busco != "off") {
        log.warn "[ASCC WARN]: run_autofilter_assembly is off, but run_btk_busco != off"
        log.warn "This will stop blobtoolkit from running unless you restart with:"
        log.warn "    `--btk_busco_run_mode mandatory`"
    }

    // Noticed a race condition, this should fix that.
    //
    run_btk_conditional.run_btk
        .map { meta, file, data -> [meta.id, meta, file] }
        .join(
            ch_autofilt_alarm_file
                .map { meta, file ->
                    [meta.id, meta, file]
                }
        )
        .map { id, ref_meta, ref_file, alarm_meta, alarm_file ->
            def merged_meta = ref_meta + alarm_meta
            [merged_meta, ref_file, alarm_file]
        }
        .set { combined_ch }


    //
    // MODULE: THIS MODULE FORMATS THE INPUT DATA IN A SPECIFIC CSV FORMAT FOR
    //          USE IN THE BTK PIPELINE
    //          EXEC MODULE PRODUCES NO VERSIONS
    //
    combined_ch
        .combine( reads_path.collect() )
        .map { meta, ref, alarm, path_list ->
            [[id:meta.id], path_list]
        }
        .set { ch_meta_reads }

    BLOBTOOLKIT_GENERATECSV (
        ch_meta_reads,
        [[],[]],
        [[],[],[]]
    )
    ch_versions     = ch_versions.mix(BLOBTOOLKIT_GENERATECSV.out.versions)


    //
    // LOGIC: STRIP THE META DATA DOWN TO id AND COMBINE ON THAT.
    //
    btk_samplesheet = BLOBTOOLKIT_GENERATECSV.out.csv
        .map{ meta, csv ->
            [[id: meta.id], csv]
        }


    //
    // So autofilter needs to be in a "Shreodingers cat" situation
    // It can either exist or not but both need to be able to run.
    // WITH AUTOFILTER
    // we can bind this file into the required inputs
    // this is to avoid a possible race condition which a generic fcs_gx (no meta)
    // will trigger btk to start running however if the PRIMARY passes AUTOFILTER
    // but HAPLO completes the other required steps
    // then HAPLO will be triggered for BTK not PRIMARY which would be correct
    // WITHOUT AUTOFILTER
    // an empty tuple [[id: "NA"], file]
    combined_input = run_btk_conditional.run_btk
        .map{ meta, file, _data ->
            [[id: meta.id], file]
        }
        .combine(btk_samplesheet, by: 0)

    combined_input
        .map{ meta, ref, samplesheet ->
            log.info("[ASCC INFO]: BTK will run for $meta\n\t| REF: ${ref}\n\t| SST: ${samplesheet}\n")
        }

    //
    // PIPELINE: PREPARE THE DATA FOR USE IN THE SANGER-TOL/BLOBTOOLKIT PIPELINE
    //              WE ARE USING THE PIPELINE HERE AS A MODULE THIS REQUIRES IT
    //              TO BE USED AS A AN INTERACTIVE JOB ON WHAT EVER EXECUTOR YOU ARE USING.
    //              This will also eventually check for the above run_btk boolean from
    //              autofilter
    SANGER_TOL_BTK (
        combined_input,
        diamond_uniprot_db_path.first(),
        nt_database_path.first(),
        diamond_uniprot_db_path.first(),
        ncbi_taxonomy_path.first(),
        reads_path.collect(),
        file("${projectDir}/assets/btk_config_files/btk_pipeline.config"),
        file("${projectDir}/assets/btk_config_files/btk_trace.config"),
        btk_lineages_path.first(),
        btk_lineages.first(),
        taxid.first(),
    )
    ch_versions     = ch_versions.mix(SANGER_TOL_BTK.out.versions)


    //-------------------------------------------------------------------------
    if (
            ( params.run_merge_datasets in run_conditionals ) &&
            ( params.run_btk_busco in run_conditionals )
        ) {
            //
            // MODULE: MERGE THE TWO BTK FORMATTED DATASETS INTO ONE DATASET FOR EASIER USE
            //
            merged_channel = ch_create_btk_dataset
                .map { meta, file -> [meta.id, [meta, file]] }
                .join(
                    SANGER_TOL_BTK.out.dataset
                        .map { meta, file ->
                            [meta.id, [meta, file]]
                    })
                .map { _id, ref_meta, ref_file, _btk_meta, btk_file ->
                    [ref_meta, ref_file, btk_file]
                }

            MERGE_BTK_DATASETS (
                merged_channel
            )
            ch_versions             = ch_versions.mix(MERGE_BTK_DATASETS.out.versions)
            busco_merge_btk         = MERGE_BTK_DATASETS.out.busco_summary_tsv
                                        .map{ meta, _tsv -> [[id: meta.id], _tsv] }
            merged_ds               = MERGE_BTK_DATASETS.out.merged_datasets
        } else {
            busco_merge_btk         = channel.of( [[:],[]] )
            merged_ds               = channel.of( [[:],[]] )

        }


    //-------------------------------------------------------------------------
    //
    // LOGIC: EACH SUBWORKFLOW OUTPUTS EITHER AN EMPTY CHANNEL OR A FILE CHANNEL DEPENDING ON THE RUN RULES
    //          SO THE RULES FOR THIS ONLY NEED TO BE A SIMPLE "DO YOU WANT IT OR NOT"
    //
    if (
        ( params.run_essentials in run_conditionals ) &&
        ( params.run_merge_datasets in run_conditionals )
    ) {

        //
        // LOGIC: JOIN CHANNELS INTO ONE BASED ON META.ID WHILST RETAINING EMPTY CHANNELS
        //
        ej_reference_tuple
            .map{meta, file -> [[id: meta.id], file]}
            .join(ej_gc_coverage
                .map{meta, file -> [[id: meta.id], file]},   remainder: true)
            .join(ch_coverage,      remainder: true)
            .join(ch_tiara,         remainder: true)
            .join(ch_kraken3,       remainder: true)
            .join(ch_blast_lineage, remainder: true)
            .join(ch_kmers,         remainder: true)
            .join(nr_hits,          remainder: true)
            .join(nr_hits,          remainder: true)
            .join(ch_create_summary,remainder: true)
            .join(busco_merge_btk,  remainder: true)
            .join(ch_fcsgx,         remainder: true)
            .filter { items ->
                def meta = items[0]
                meta != null &&
                meta != [] &&
                !(meta instanceof Map && (meta.id == null || meta.isEmpty()))
            }
            .map { items ->
                // Replace null values with placeholder file
                items.withIndex().collect { item, index ->
                    if (item == null) {
                        getEmptyPlaceholder(index)
                    } else if (item instanceof List && item.isEmpty()) {
                        getEmptyPlaceholder(index)
                    } else {
                        item
                    }
                }
            }
            .set{ merge_input_channel}

        ASCC_MERGE_TABLES (
            merge_input_channel
        )
        ch_versions             = ch_versions.mix(ASCC_MERGE_TABLES.out.versions)

        merged_table            = ASCC_MERGE_TABLES.out.merged_table
                                    .map{ meta, _file -> [[id: meta.id ], _file] }

        merged_extended_table   = ASCC_MERGE_TABLES.out.extended_table
        merged_phylum_count     = ASCC_MERGE_TABLES.out.phylum_counts
                                    .map{ meta, _file -> [[id: meta.id], _file] }
    } else {
        merged_table            = channel.of( [[:],[]] )
        merged_extended_table   = channel.empty()
        merged_phylum_count     = channel.of( [[:],[]] )
    }


    //-------------------------------------------------------------------------
    //
    // SUBWORKFLOW: GENERATE HTML REPORT (minimal wiring, opt-in)
    //              Gate with params.run_html_report to avoid altering default behavior.
    //

    // Params file
    ch_params_file      = params.params_file ? channel.fromPath(params.params_file) : channel.value([])

    // Templates and CSS
    ch_jinja_templates = channel.fromPath("${baseDir}/assets/templates/*.jinja").collect()
    ch_css_files       = channel.fromPath("${baseDir}/assets/css/*.css").collect()

    GENERATE_HTML_REPORT_WORKFLOW (
        ch_barcode_check,
        ch_fcsadapt,
        ej_trailing_ns,
        ch_vecscreen,
        ch_autofilt_fcs_tiara,
        merged_table,
        merged_phylum_count,
        ch_kmers_results,
        ej_reference_tuple.filter{ meta, file ->
            params.run_html_report in run_conditionals
        },
        ej_fasta_sanitation_log,
        ej_fasta_filter_log,
        ch_params_file,
        ch_fcsgx_report,
        ch_fcsgx_taxonomy,
        ch_create_btk_dataset
    )
    ch_versions             = ch_versions.mix(GENERATE_HTML_REPORT_WORKFLOW.out.versions)

    emit:
    essential_reference         = ej_reference_tuple
    essential_genome_file       = ej_dot_genome
    essential_gc_cov            = ej_gc_coverage

    kmer_data                   = ch_kmers

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

    organellar_blast_mito       = ch_mito
    organellar_blast_chloro     = ch_chloro

    pacbio_barcode_files        = ch_barcode_check // This is a collection of (params.barcode * [meta, file])

    ascc_merged_table           = merged_table
    ascc_merged_table_extended  = merged_extended_table
    ascc_merged_table_phylum_c  = merged_phylum_count

    merged_btk_ds_datasets      = merged_ds
    merged_btk_ds_busco_summary = busco_merge_btk

    // THESE ONES DON'T RELY ON THE NORMAL IF ELSE STRUCTURE OF THE OTHER
    // SUBWORKFLOWS SO THERE'S NO "BACKUP" CHANNEL.
    // sanger_tol_btk_dataset      = SANGER_TOL_BTK.out.dataset
    // sanger_tol_btk_plots        = SANGER_TOL_BTK.out.plots
    // sanger_tol_btk_summary_json = SANGER_TOL_BTK.out.summary_json
    // sanger_tol_btk_busco_data   = SANGER_TOL_BTK.out.busco_data
    // sanger_tol_btk_multiqc      = SANGER_TOL_BTK.out.multiqc_report
    // sanger_tol_btk_pipeline_info= SANGER_TOL_BTK.out.pipeline_info

    // generate_samplesheet_csv    = GENERATE_SAMPLESHEET.out.csv

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
