/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CREATE_BTK_DATASET                            } from '../modules/local/blobtoolkit/create_dataset/main'
include { AUTOFILTER_AND_CHECK_ASSEMBLY                 } from '../modules/local/autofilter/autofilter/main'

include { ESSENTIAL_JOBS                                } from '../subworkflows/local/essential_jobs/main'
include { EXTRACT_TIARA_HITS                            } from '../subworkflows/local/extract_tiara_hits/main'
include { EXTRACT_NT_BLAST                              } from '../subworkflows/local/extract_nt_blast/main'
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

workflow ASCC_ORGANELLAR {

    take:
    ch_samplesheet          // channel: samplesheet read in from --input
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
    reads_type

    main:
    ch_versions = Channel.empty()


    //
    // LOGIC: PRETTY NOTIFICATION OF FILES AT STAGE
    //
    ch_samplesheet
        .map { meta, sample ->
            log.info "[ASCC info] ORGANELLAR WORKFLOW:\n\t-- $meta\n\t-- $sample\n"
        }


    //
    // SUBWORKFLOW: RUNS FILTER_FASTA, GENERATE .GENOME, CALCS GC_CONTENT AND FINDS RUNS OF N's
    //
    if ( params.run_essentials == "both" || params.run_essentials == "organellar" ) {
        ESSENTIAL_JOBS(
            ch_samplesheet
        )
        ch_versions             = ch_versions.mix(ESSENTIAL_JOBS.out.versions)
        reference_tuple_from_GG = ESSENTIAL_JOBS.out.reference_tuple_from_GG
        reference_tuple_w_seqkt = ESSENTIAL_JOBS.out.reference_with_seqkit
        ej_dot_genome           = ESSENTIAL_JOBS.out.dot_genome
        ej_gc_coverage          = ESSENTIAL_JOBS.out.gc_content_txt

    } else {
        log.warn("[ASCC warn] MAKE SURE YOU ARE AWARE YOU ARE SKIPPING ESSENTIAL JOBS, THIS INCLUDES BREAKING SCAFFOLDS OVER 1.9GB, FILTERING N\'s AND GC CONTENT REPORT (THIS WILL BREAK OTHER PROCESSES AND SHOULD ONLY BE RUN WITH `--include essentials`)")

        reference_tuple_from_GG = ch_samplesheet
        ej_dot_genome           = Channel.empty()
        ej_gc_coverage          = Channel.empty()
        reference_tuple_w_seqkt = Channel.empty()
    }


    //
    // SUBWORKFLOW: EXTRACT RESULTS HITS FROM TIARA
    //
    if ( params.run_tiara == "both" || params.run_tiara == "organellar" ) {
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
    // SUBWORKFLOW: IDENTITY PACBIO BARCODES IN INPUT DATA
    //
    if ( params.run_pacbio_barcodes == "both" || params.run_pacbio_barcodes == "organellar" ) {

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
    if ( params.run_fcs_adaptor == "both" || params.run_fcs_adaptor == "organellar" ) {
        RUN_FCSADAPTOR (
            reference_tuple_from_GG
        )

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

    if ( (params.run_fcsgx == "both" || params.run_fcsgx == "organellar") & !params.fcs_override) {

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
        ch_fcsgx            = RUN_FCSGX.out.fcsgxresult
                                .map { it ->
                                    [[id: it[0].id, process: "FCSGX result"], it[1]]
                                }
                                .ifEmpty { [[],[]] }

    } else if (params.fcs_override) {
        log.info("[ASCC info] Overriding Internal FCSGX")
        ch_fcsgx         = fcs_ss
        ch_fcsgx.view{"[ASCC info] OVERRIDDEN_FCSGX: $it"}
    } else {
        ch_fcsgx         = Channel.of( [[],[]] )
    }


    //
    // SUBWORKFLOW: CALCULATE AVERAGE READ COVERAGE
    //
    if ( params.run_coverage == "both" || params.run_coverage == "genomic" ) {

        RUN_READ_COVERAGE (
            reference_tuple_from_GG, // Again should this be the validated fasta?
            reads,
            reads_type,
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
    if ( params.run_vecscreen == "both" || params.run_vecscreen == "genomic" ) {

        RUN_VECSCREEN (
            reference_tuple_from_GG, // Again should this be the validated fasta?
            vecscreen_database_path.first()
        )
        ch_versions         = ch_versions.mix(RUN_VECSCREEN.out.versions)
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


    //
    // LOGIC: WE NEED TO MAKE SURE THAT THE INPUT SEQUENCE IS OF AT LEAST LENGTH OF params.seqkit_window
    //
    valid_length_fasta = reference_tuple_w_seqkt
        //
        // NOTE: Here we are using the un-filtered genome, any filtering may (accidently) cause an empty fasta
        //
        .map{ meta, file ->
            tuple(
                [
                    id: meta.id,
                    sliding: meta.sliding,
                    window: meta.window,
                    seq_count: CountFastaLength(file)
                ],
                file
            )
        }
        .filter { meta, file ->
                    meta.seq_count >= params.seqkit_window
        }

    valid_length_fasta
        .map{ meta, file ->
            log.info "[ASCC info] Running BLAST (NT, DIAMOND, NR) on VALID ORGANELLE: \n\t-- ${meta.id}'s sequence ($meta.seq_count bases) is >= seqkit_window $params.seqkit_window\n"
        }

    //
    // LOGIC: THIS CONDITIONAL SHOULD EXECUTE THE PROCESS WHEN:
    //          INCLUDE STEPS ARE EITHER nt_blast AND all
    //              _AS WELL AS_
    //          EXCLUDE _NOT_ CONTAINING nt_blast AND THE valid_length_fasta IS NOT EMPTY
    //
    if ( params.run_nt_blast == "both" || params.run_nt_blast == "genomic" ) {

        //
        //SUBWORKFLOW: EXTRACT RESULTS HITS FROM NT-BLAST
        //

        EXTRACT_NT_BLAST (
            valid_length_fasta,
            nt_database_path.first(),
            ncbi_ranked_lineage_path.first()
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
    if ( params.run_nr_diamond == "both" || params.run_nr_diamond == "genomic" ) {

        NR_DIAMOND (
            valid_length_fasta,
            diamond_nr_db_path.first()
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
    //
    //qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms sphylums salltitles
    if ( params.run_uniprot_diamond == "both" || params.run_uniprot_diamond == "genomic" ) {

        UP_DIAMOND (
            valid_length_fasta,
            diamond_uniprot_db_path.first()
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


    if ( params.run_create_btk_dataset == "both" || params.run_create_btk_dataset == "genomic" ) {

        //
        // LOGIC: FOUND RACE CONDITION EFFECTING LONG RUNNING JOBS
        //          AND INPUT TO HERE ARE NOW MERGED AND MAPPED
        //          EMPTY CHANNELS ARE CHECKED AND DEFAULTED TO [[],[]]
        //
        ch_organellar_cbtk_input = reference_tuple_from_GG
            .map{ it -> tuple([
                id: it[0].id,
                taxid: it[0].taxid,
                sci_name: it[0].sci_name,
                process: "REFERENCE"], it[1])
            }
            .mix(
                ej_dot_genome.map{ it -> tuple([id: it[0].id, process: "GENOME"], it[1])},
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

        ch_autofilt_alarm_file  = AUTOFILTER_AND_CHECK_ASSEMBLY.out.alarm_file
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
        ch_autofilt_alarm_file  = Channel.empty()
        ch_autofilt_removed_seqs= Channel.empty()
        ch_autofilt_assem       = Channel.empty()
        ch_autofilt_indicator   = Channel.empty()
        ch_autofilt_fcs_tiara   = Channel.empty()
        ch_autofilt_raw_report  = Channel.empty()
    }

    emit:

    essential_reference         = reference_tuple_from_GG
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

    create_btk_ds_dataset       = create_btk_dataset
    create_btk_ds_create_smry   = create_summary

    kraken2_classified          = ch_kraken1
    kraken2_report              = ch_kraken2
    kraken2_lineage             = ch_kraken3

    vecscreen_contam            = ch_vecscreen

    tiara_output                = ch_tiara

    versions                    = ch_versions

}

//
// Function: this is to count the length of ONLY the fasta sequence
//
// @param input_file: path
// @return int
def CountFastaLength(input_file) {
    int counter = 0;
    def list_lines = new File(input_file.toString()).text.readLines()
    for (i in list_lines) {
        if (i[0] != ">") {
            counter = counter + i.length()
        }
    }
    return counter;
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
