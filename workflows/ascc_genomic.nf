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
    fcs_ov                  // params.fcs_override
    fcs_ss                  //
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
            log.info "[ASCC info] GENOMIC WORKFLOW:\n\t-- $meta\n\t-- $sample\n"
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

    } else {
        log.warn("[ASCC warn] MAKE SURE YOU ARE AWARE YOU ARE SKIPPING ESSENTIAL JOBS, THIS INCLUDES BREAKING SCAFFOLDS OVER 1.9GB, FILTERING N\'s AND GC CONTENT REPORT (THIS WILL BREAK OTHER PROCESSES AND SHOULD ONLY BE RUN WITH `--include essentials`)")

        reference_tuple_from_GG = ch_samplesheet
        ej_dot_genome           = Channel.empty()
        ej_gc_coverage          = Channel.empty()
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
        ch_barcode_check    = PACBIO_BARCODE_CHECK.out.filtered.collect()
        ch_versions         = ch_versions.mix(PACBIO_BARCODE_CHECK.out.versions)

    } else {
        ch_barcode_check    = Channel.empty()
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
    if ( (params.run_fcsgx == "both" || params.run_fcsgx == "genomic") && !params.fcs_override ) {

        joint_channel = reference_tuple_from_GG
            .combine(fcs_db)
            .combine(taxid)
            .combine(ncbi_ranked_lineage_path)
            .multiMap { meta, ref, db, tax_id, tax_path ->
                meta = [id: meta.id, taxid: meta.taxid]
                reference: [meta, ref]
                fcs_db_path: db
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
                                .map { meta, file ->
                                    [[id: meta.id, process: "FCSGX result"], file]
                                }
                                .ifEmpty { [[],[]] }

    } else if ( params.fcs_override ) {
        log.info("[ASCC info] Overriding Internal FCSGX")
        ch_fcsgx            = fcs_ss

        ch_fcsgx.view{"OVERRIDDEN_FCSGX: $it"}
    } else {
        ch_fcsgx            = Channel.of( [[],[]] )
    }


    //
    // SUBWORKFLOW: CALCULATE AVERAGE READ COVERAGE
    //
    //

    if ( params.run_coverage == "both" || params.run_coverage == "genomic" ) {
        RUN_READ_COVERAGE (
            reference_tuple_from_GG,
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
                                .ifEmpty { [[],[]] }

        ch_bam              = RUN_READ_COVERAGE.out.bam_ch
                                .map { meta, file ->
                                    tuple([id: meta.id, process: "Mapped Bam"], file)
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
            scientific_name

        )
        ch_versions             = ch_versions.mix(CREATE_BTK_DATASET.out.versions)

        create_summary          = CREATE_BTK_DATASET.out.create_summary.map{ it -> tuple([id: it[0].id, process: "C_BTK_SUM"], it[1])}
        create_btk_dataset      = CREATE_BTK_DATASET.out.btk_datasets
    } else {
        create_summary          = Channel.empty()
        create_btk_dataset      = Channel.empty()
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
                    def new_meta = [id: meta.id, taxid: thetaxid]
                    reference:  tuple(new_meta, ref)
                    tiara_file: tuple(new_meta, tiara)
                    fcs_file:   tuple(new_meta, fcs)
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
            .map { meta, file -> [meta, file.text.trim()] }
            .branch { meta, data ->
                log.info("[ASCC info] Run for ${meta.id} has:\n\t${data}\n")

                run_btk     : data.contains("YES_ABNORMAL_CONTAMINATION") ? tuple(meta, "YES") : Channel.empty()
                dont_run    : true
            }

        btk_bool_run_btk        = btk_bool.run_btk

        ch_autofilt_alarm_file   = AUTOFILTER_AND_CHECK_ASSEMBLY.out.alarm_file
            .map{ meta, file ->
                tuple(
                    [id: meta.id], file
                )
            }

        ch_autofilt_fcs_tiara   = AUTOFILTER_AND_CHECK_ASSEMBLY.out.fcs_tiara_summary
        ch_autofilt_removed_seqs= AUTOFILTER_AND_CHECK_ASSEMBLY.out.removed_seqs
        ch_autofilt_raw_report  = AUTOFILTER_AND_CHECK_ASSEMBLY.out.raw_report

        ch_versions             = ch_versions.mix(AUTOFILTER_AND_CHECK_ASSEMBLY.out.versions)
    } else {
        btk_bool_run_btk        = Channel.of([[id: "NA"], "false"])
        ch_autofilt_alarm_file  = Channel.empty()
        ch_autofilt_removed_seqs= Channel.empty()
        ch_autofilt_assem       = Channel.empty()
        ch_autofilt_indicator   = Channel.empty()
        ch_autofilt_fcs_tiara   = Channel.empty()
        ch_autofilt_raw_report  = Channel.empty()
    }


    //
    // LOGIC: DETERMINE WHETHER BLOBTOOLKIT SHOULD RUN BASED ON CONDITIONALS
    //         - ALWAYS RUN IF params.btk_busco_run_mode == "mandatory" AND BTK

    run_btk_conditional = reference_tuple_from_GG
        | map { meta, file ->
                tuple([id: meta.id, taxid: meta.taxid] , file)
            }
        // below is combined into the tuple to enforce the block to only run when channel is present.
        | combine ( btk_bool_run_btk
                        .map{ meta, data ->
                            def joined_content = data
                                    .replaceAll(/\s*\|\s*/, "-")     // Replace " | " with "-"
                                    .replaceAll(/\s+/, "-")          // Replace remaining spaces with "-"
                                    .replaceAll(/_+/, "_")           // Keep underscores as they are
                                    .replaceAll(/-+/, "-")           // Clean up multiple dashes
                            tuple(
                                [id: meta.id, taxid: meta.taxid],
                                joined_content
                            )
                        },
                by: [0]
            )
        | branch { meta, assembly, data ->
            def btk_requested           = params.run_btk_busco == "both" || params.run_btk_busco == "genomic"
            def autofilter_requested    = params.run_autofilter_assembly == "both" || params.run_autofilter_assembly == "genomic"

            def ignore_autofilter       = params.btk_busco_run_mode == "mandatory" && btk_requested
            def not_mandatory_btk       = params.btk_busco_run_mode == "conditional" && autofilter_requested && btk_requested && data.contains("YES")

            run_btk: (ignore_autofilter || not_mandatory_btk)
            skip_btk: true
        }

    run_btk_conditional.skip_btk
        .map { meta, file, data ->
            log.warn "[ASCC WARNING]: CONTAMINATION THRESHOLD NOT MET"
            log.warn "\t- SKIPPING BLOBTOOLKIT FOR: $meta.id"
            log.warn "\t- You can verify here: $file"
        }

    if (params.run_autofilter_assembly == "off" && params.run_btk_busco != "off") {
        log.info "[ASCC info] run_autofilter_assembly is off, but run_btk_busco != off"
        log.info "This will stop blobtoolkit from running unless you restart with:"
        log.info "    `--btk_busco_run_mode mandatory`"
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
    //
    //
    GENERATE_SAMPLESHEET (
        combined_ch,
        reads_path.first(),
        reads_layout.first()
    )
    ch_versions         = ch_versions.mix(GENERATE_SAMPLESHEET.out.versions)

    //
    // LOGIC: STRIP THE META DATA DOWN TO id AND COMBINE ON THAT.
    //
    btk_samplesheet = GENERATE_SAMPLESHEET.out.csv
        .map{ meta, csv, indicator ->
            tuple(
                [ id: meta.id ],
                csv,
                indicator
            )
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
        .map{
            it -> tuple(
                [id:it[0].id],
                it[1]
            )
        }
        .combine(btk_samplesheet, by: 0)

    combined_input
        .map{ meta, ref, samplesheet, alarms ->
            log.info("[ASCC info] BTK will run for $meta\n\t| REF: ${ref}\n\t| SST: ${samplesheet}\n\t| ALM: ${alarms}\n")
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
        reads_path.first(),
        file("${projectDir}/assets/btk_config_files/btk_pipeline.config"),
        btk_lineages_path.first(),
        btk_lineages.first(),
        taxid.first(),
    )
    ch_versions             = ch_versions.mix(SANGER_TOL_BTK.out.versions)


if (
        ( params.run_merge_datasets == "both" || params.run_merge_datasets == "genomic" ) &&
        ( params.run_btk_busco == "both" || params.run_btk_busco == "genomic" )
    ) {
        //
        // MODULE: MERGE THE TWO BTK FORMATTED DATASETS INTO ONE DATASET FOR EASIER USE
        //
        merged_channel = create_btk_dataset
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
        merged_ds               = MERGE_BTK_DATASETS.out.merged_datasets
    } else {
        busco_merge_btk         = Channel.empty()
        merged_ds               = Channel.empty()

    }


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

        merged_table            = ASCC_MERGE_TABLES.out.merged_table
        merged_extended_table   = ASCC_MERGE_TABLES.out.extended_table
        merged_phylum_count     = ASCC_MERGE_TABLES.out.phylum_counts

    } else {
        merged_table            = Channel.empty()
        merged_extended_table   = Channel.empty()
        merged_phylum_count     = Channel.empty()
    }

    emit:
    essential_reference         = reference_tuple_from_GG
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

    create_btk_ds_dataset       = create_btk_dataset
    create_btk_ds_create_smry   = create_summary

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
