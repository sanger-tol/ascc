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
include { GENERATE_SAMPLESHEET                          } from '../modules/local/blobtoolkit/generate_samplesheet/main'
include { NEXTFLOW_RUN as SANGER_TOL_BTK_CASCADE        } from '../modules/local/run/main'

include { ESSENTIAL_JOBS                                } from '../subworkflows/local/essential_jobs/main'
include { GET_KMERS_PROFILE                             } from '../subworkflows/local/get_kmers_profile/main'
include { EXTRACT_TIARA_HITS                            } from '../subworkflows/local/extract_tiara_hits/main'
include { EXTRACT_NT_BLAST                              } from '../subworkflows/local/extract_nt_blast/main'
include { ORGANELLAR_BLAST as PLASTID_ORGANELLAR_BLAST  } from '../subworkflows/local/organellar_blast/main'
include { ORGANELLAR_BLAST as MITO_ORGANELLAR_BLAST     } from '../subworkflows/local/organellar_blast/main'
include { PACBIO_BARCODE_CHECK                          } from '../subworkflows/local/pacbio_barcode_check/main'
include { TRAILINGNS_CHECK                              } from '../subworkflows/local/trailingns_check/main'
include { RUN_READ_COVERAGE                             } from '../subworkflows/local/run_read_coverage/main'
include { RUN_VECSCREEN                                 } from '../subworkflows/local/run_vecscreen/main'
include { RUN_NT_KRAKEN                                 } from '../subworkflows/local/run_nt_kraken/main'
include { RUN_FCSGX                                     } from '../subworkflows/local/run_fcsgx/main'
include { RUN_FCSADAPTOR                                } from '../subworkflows/local/run_fcsadaptor/main'
include { RUN_DIAMOND as NR_DIAMOND                     } from '../subworkflows/local/run_diamond/main'
include { RUN_DIAMOND as UP_DIAMOND                     } from '../subworkflows/local/run_diamond/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASCC_GENOMIC {

    take:
    ch_samplesheet          // channel: samplesheet read in from --input
    organellar_genomes      // channel: tuple(meta, reference)
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

    main:
    ch_versions = Channel.empty()

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

    if ( params.run_essentials == "both" || params.run_essentials == "genomic" ) {
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

        reference_tuple_from_GG = ch_samplesheet
        ej_dot_genome           = Channel.empty()
        ej_gc_coverage          = Channel.empty()
        reference_tuple_w_seqkt = Channel.empty()
    }


    if ( params.run_kmers == "both" || params.run_kmers == "genomic" ) {
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

        //
        // LOGIC: AT THIS POINT THE META CONTAINS JUNK THAT CAN 'CONTAMINATE' MATCHES,
        //          SO STRIP IT DOWN AND ADD PROCESS_NAME BEFORE USE
        //
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
    if ( params.run_tiara == "both" || params.run_tiara == "genomic" ) {
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
    if ( params.run_nt_blast == "both" || params.run_nt_blast == "genomic" ) {
        //
        // NOTE: ch_nt_blast needs to be set in two places incase it
        //          fails during the run
        //
        EXTRACT_NT_BLAST (
            reference_tuple_from_GG,
            nt_database_path.first(),
            ncbi_ranked_lineage_path.first()
        )
        ch_versions         = ch_versions.mix(EXTRACT_NT_BLAST.out.versions)

        //
        // LOGIC: AT THIS POINT THE META CONTAINS JUNK THAT CAN 'CONTAMINATE' MATCHES,
        //          SO STRIP IT DOWN AND ADD PROCESS_NAME BEFORE USE
        //
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
    if ( params.run_nr_diamond == "both" || params.run_nr_diamond == "genomic" ) {

        NR_DIAMOND (
            reference_tuple_from_GG,
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
    if ( params.run_uniprot_diamond == "both" || params.run_uniprot_diamond == "genomic" ) {

        UP_DIAMOND (
            reference_tuple_from_GG,
            diamond_uniprot_db_path.first()
        )
        ch_versions         = ch_versions.mix(UP_DIAMOND.out.versions)

        //
        // LOGIC: AT THIS POINT THE META CONTAINS JUNK THAT CAN 'CONTAMINATE' MATCHES,
        //          SO STRIP IT DOWN AND ADD PROCESS_NAME BEFORE USE
        //
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


    if ( params.run_organellar_blast == "both" || params.run_organellar_blast == "genomic" ) {
        //
        // LOGIC: CHECK WHETHER THERE IS A MITO AND BRANCH
        //
        organellar_check = organellar_genomes
            .branch { meta, assembly ->
                mito:       meta.assembly_type == "MITO"
                plastid:    meta.assembly_type == "PLASTID"
                invalid:    true    // if value but not of the above conditions
            }


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


        //
        // LOGIC: AT THIS POINT THE META CONTAINS JUNK THAT CAN 'CONTAMINATE' MATCHES,
        //          SO STRIP IT DOWN AND ADD PROCESS_NAME BEFORE USE
        //
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
    if ( params.run_pacbio_barcodes == "both" || params.run_pacbio_barcodes == "genomic" ) {

        reference_tuple_from_GG
            .combine(pacbio_database)
            .multiMap{
                ref_meta, ref_data, pdb_meta, pdb_data ->
                    reference: [ref_meta, ref_data]
                    pacbio_db: [pdb_meta, pdb_data]
            }
            .set { duplicated_db }

        PACBIO_BARCODE_CHECK (
            duplicated_db.reference,
            params.pacbio_barcode_names,
            duplicated_db.pacbio_db
        )

        ch_versions         = ch_versions.mix(PACBIO_BARCODE_CHECK.out.versions)
    }


    //
    // SUBWORKFLOW: RUN FCS-ADAPTOR TO IDENTIDY ADAPTOR AND VECTORR CONTAMINATION
    //
    if ( params.run_fcs_adaptor == "both" || params.run_fcs_adaptor == "genomic" ) {
        RUN_FCSADAPTOR (
            reference_tuple_from_GG
        )
        ch_versions         = ch_versions.mix(RUN_FCSADAPTOR.out.versions)

        //
        // LOGIC: AT THIS POINT THE META CONTAINS JUNK THAT CAN 'CONTAMINATE' MATCHES,
        //          SO STRIP IT DOWN BEFORE USE, WE ALSO MERGE THE OUTPUT TOGETHER FOR SIMPLICITY
        //
        ch_fcsadapt = RUN_FCSADAPTOR.out.ch_euk
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

    } else {
        ch_fcsadapt         = Channel.empty()
    }


    //
    // SUBWORKFLOW: RUN FCS-GX TO IDENTIFY CONTAMINATION IN THE ASSEMBLY
    //
    if ( params.run_fcsgx == "both" || params.run_fcsgx == "genomic" ) {

        joint_channel = reference_tuple_from_GG
            .combine(fcs_db)
            .combine(taxid)
            .combine(ncbi_ranked_lineage_path)
            .multiMap { meta, ref, db, tax_id, tax_path ->
                reference: [meta, tax_id, ref]
                fcs_db_path: db
                taxid_val: tax_id
                ncbi_tax_path: tax_path
            }

        RUN_FCSGX (
            joint_channel.reference,
            joint_channel.fcs_db_path,
            joint_channel.ncbi_tax_path
        )
        ch_versions         = ch_versions.mix(RUN_FCSGX.out.versions)

        //
        // LOGIC: AT THIS POINT THE META CONTAINS JUNK THAT CAN 'CONTAMINATE' MATCHES,
        //          SO STRIP IT DOWN AND ADD PROCESS_NAME BEFORE USE
        //
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
    if ( params.run_coverage == "both" || params.run_coverage == "genomic" ) {
        RUN_READ_COVERAGE (
            reference_tuple_from_GG,
            reads.first(),
            reads_type.first(),
        )
        ch_versions         = ch_versions.mix(RUN_READ_COVERAGE.out.versions)

        //
        // LOGIC: AT THIS POINT THE META CONTAINS JUNK THAT CAN 'CONTAMINATE' MATCHES,
        //          SO STRIP IT DOWN AND ADD PROCESS_NAME BEFORE USE
        //
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
    if ( params.run_vecscreen == "both" || params.run_vecscreen == "genomic" ) {
        RUN_VECSCREEN (
            reference_tuple_from_GG,
            vecscreen_database_path.first()
        )
        ch_versions         = ch_versions.mix(RUN_VECSCREEN.out.versions)

        //
        // LOGIC: AT THIS POINT THE META CONTAINS JUNK THAT CAN 'CONTAMINATE' MATCHES,
        //          SO STRIP IT DOWN AND ADD PROCESS_NAME BEFORE USE
        //
        ch_vecscreen        = RUN_VECSCREEN.out.vecscreen_contam
                                .map { it ->
                                    [[id: it[0].id, process: "Vecscreen"], it[1]]
                                }
                                .ifEmpty { [[],[]] }
    } else {
        ch_vecscreen        = Channel.empty()
    }


    //
    // SUBWORKFLOW: RUN THE KRAKEN CLASSIFIER
    //
    if ( params.run_kraken == "both" || params.run_kraken == "genomic" ) {

        RUN_NT_KRAKEN(
            reference_tuple_from_GG,
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
        ch_kraken1 = Channel.empty()
        ch_kraken2 = Channel.empty()
        ch_kraken3 = Channel.empty()

    }


    if ( params.run_create_btk_dataset == "both" || params.run_create_btk_dataset == "genomic" ) {

        //
        // LOGIC: FOUND RACE CONDITION EFFECTING LONG RUNNING JOBS
        //          AND INPUT TO HERE ARE NOW MERGED AND MAPPED
        //          EMPTY CHANNELS ARE CHECKED AND DEFAULTED TO [[],[]]
        //
        ch_genomic_cbtk_input = reference_tuple_from_GG
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
            ncbi_taxonomy_path.first(),
            scientific_name.first()

        )
        ch_versions             = ch_versions.mix(CREATE_BTK_DATASET.out.versions)

        create_summary          = CREATE_BTK_DATASET.out.create_summary.map{ it -> tuple([id: it[0].id, process: "C_BTK_SUM"], it[1])}
    } else {
        create_summary          = Channel.empty()
    }


    //
    // LOGIC: AUTOFILTER ASSEMBLY BY TIARA AND FCSGX RESULTS SO THE SUBWORKLOW CAN EITHER BE TRIGGERED BY THE VALUES tiara, fcs-gx, autofilter_assemlby AND EXCLUDE STEPS NOT CONTAINING autofilter_assembly
    //          OR BY include_steps CONTAINING ALL AND EXCLUDE NOT CONTAINING autofilter_assembly.
    //
    if (
        ( params.run_tiara == "both" || params.run_tiara == "genomic" ) &&
        ( params.run_fcsgx == "both" || params.run_fcsgx == "genomic" ) &&
        ( params.run_autofilter_assembly == "both" || params.run_autofilter_assembly == "genomic" )
    ) {
        //
        // LOGIC: FILTER THE INPUT FOR THE AUTOFILTER STEP
        //          - We can't just combine on meta.id as some of the Channels have other data
        //              in there too so we just sanitise, and _then_ combine on 0, and
        //              _then_ add back in the taxid as we need that for this process.
        //              Thankfully taxid is a param so easy enough to add back in.
        //                  Actually, it just makes more sense to passs in as its own channel.
        //

        autofilter_input_formatted = reference_tuple_from_GG
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
                ncbi_ranked_lineage_path
            )
            .combine(
                taxid
            )
            .multiMap{
                meta, ref, tiara, fcs, ncbi, thetaxid ->
                    reference:  tuple([id: meta.id, taxid: thetaxid], ref)
                    tiara_file: tuple(meta, tiara)
                    fcs_file:   tuple(meta, fcs)
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

        //
        // LOGIC: BRANCH THE CHANNEL ON WHETHER OR NOT THERE IS ABNORMAL CONTAMINATION IN THE
        //          OUTPUT FILE.
        //          CHANGE OUTPUT NAME TO BE REFERENCE NAME AND THEN ALARM FILE
        btk_bool = AUTOFILTER_AND_CHECK_ASSEMBLY.out.alarm_file
            .map { file -> file.text.trim() }
            .branch { it ->
                run_btk     : it.contains("YES_ABNORMAL_CONTAMINATION")
                dont_run    : true
            }

        ch_versions             = ch_versions.mix(AUTOFILTER_AND_CHECK_ASSEMBLY.out.versions)
    } else {
        ch_autofilt_assem       = Channel.empty()
        ch_autofilt_indicator   = Channel.empty()
    }


    //
    // LOGIC: DETERMINE WHETHER BLOBTOOLKIT SHOULD RUN BASED ON CONDITIONALS
    //         - ALWAYS RUN IF params.btk_busco_run_mode == "mandatory" AND BTK

    run_btk_conditional = reference_tuple_from_GG
        | branch { meta, assembly ->
            def btk_requested           = params.run_btk_busco == "both" || params.run_btk_busco == "genomic"
            def autofilter_requested    = params.run_autofilter_assembly == "both" || params.run_autofilter_assembly == "genomic"

            def ignore_autofilter       = params.btk_busco_run_mode == "mandatory" && btk_requested
            def not_mandatory_btk       = params.btk_busco_run_mode == "conditional" && autofilter_requested && btk_requested && btk_bool.run_btk

            run_btk: (ignore_autofilter || not_mandatory_btk)
            skip_btk: true
        }

    run_btk_conditional.skip_btk
        .map { meta, file ->
            log.warn "BLOBTOOLKIT: SKIPPING BLOBTOOLKIT FOR: [$meta, $file]"
        }

    if (params.run_autofilter_assembly == "off" && params.run_btk_busco != "off") {
        log.info "run_autofilter_assembly is off, but run_btk_busco != off"
        log.info "This will stop blobtoolkit from running unless you restart with:"
        log.info "    `--btk_busco_run_mode mandatory`"
    }


    //
    // MODULE: THIS MODULE FORMATS THE INPUT DATA IN A SPECIFIC CSV FORMAT FOR
    //          USE IN THE BTK PIPELINE
    //
    GENERATE_SAMPLESHEET (
        run_btk_conditional.run_btk,
        reads_path.first(),
        reads_layout.first()
    )
    ch_versions         = ch_versions.mix(GENERATE_SAMPLESHEET.out.versions)

    //
    // LOGIC: STRIP THE META DATA DOWN TO id AND COMBINE ON THAT.
    //
    btk_samplesheet = GENERATE_SAMPLESHEET.out.csv
        .map{ meta, csv ->
            tuple(
                [ id: meta.id ],
                csv
            )
        }

    combined_input = run_btk_conditional.run_btk
        .map{ it -> tuple([id:it[0].id], it[1])}
        .combine(btk_samplesheet, by: 0)
        .multiMap { meta_1, ref, csv ->
            reference: [meta_1, ref]
            samplesheet: csv
        }

    //
    // PIPELINE: PREPARE THE DATA FOR USE IN THE SANGER-TOL/BLOBTOOLKIT PIPELINE
    //              WE ARE USING THE PIPELINE HERE AS A MODULE THIS REQUIRES IT
    //              TO BE USED AS A AN INTERACTIVE JOB ON WHAT EVER EXECUTOR YOU ARE USING.
    //              This will also eventually check for the above run_btk boolean from
    //              autofilter
    SANGER_TOL_BTK (
        combined_input.reference,
        combined_input.samplesheet,
        diamond_uniprot_db_path.first(),
        nt_database_path.first(),
        diamond_uniprot_db_path.first(),
        ncbi_taxonomy_path.first(),
        reads_path.first(),
        btk_lineages_path.first(),
        btk_lineages.first(),
        taxid.first(),
    )
    ch_versions             = ch_versions.mix(SANGER_TOL_BTK.out.versions)


    //
    // MODULE: MERGE THE TWO BTK FORMATTED DATASETS INTO ONE DATASET FOR EASIER USE
    //
    merged_channel = CREATE_BTK_DATASET.out.btk_datasets
        .map { meta, file -> [meta.id, [meta, file]] }
        .join(
            SANGER_TOL_BTK.out.dataset
                .map { meta, file ->
                    [meta.id, [meta, file]]
            })
        .map { id, ref, btk -> [ref[0], ref[1], btk[1]] }

    MERGE_BTK_DATASETS (
        merged_channel
    )
    ch_versions             = ch_versions.mix(MERGE_BTK_DATASETS.out.versions)
    busco_merge_btk         = MERGE_BTK_DATASETS.out.busco_summary_tsv


    //
    // LOGIC: EACH SUBWORKFLOW OUTPUTS EITHER AN EMPTY CHANNEL OR A FILE CHANNEL DEPENDING ON THE RUN RULES
    //          SO THE RULES FOR THIS ONLY NEED TO BE A SIMPLE "DO YOU WANT IT OR NOT"
    //
    if (
        ( params.run_essentials == "both" || params.run_essentials == "genomic" ) &&
        ( params.run_merge_datasets == "both" || params.run_merge_datasets == "genomic" )
    ) {

        //
        // LOGIC: FOUND RACE CONDITION EFFECTING LONG RUNNING JOBS
        //          AND INPUT TO HERE ARE NOW MERGED AND MAPPED
        //          EMPTY CHANNELS ARE CHECKED AND DEFAULTED TO [[],[]]
        //
        ascc_merged_data = ej_gc_coverage
            .map{ it -> tuple([
                id: it[0].id,
                process: "GC_COV"], it[1])
            }
            .mix(
                ej_dot_genome.map{ it -> tuple([id: it[0].id, process: "GENOME"], it[1])},
                create_summary,
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

    emit:
    versions                    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
