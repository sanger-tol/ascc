/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ASCC GENOMIC ANALYSIS WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GET_KMERS_PROFILE                             } from '../subworkflows/local/get_kmers_profile/main'
include { EXTRACT_TIARA_HITS                            } from '../subworkflows/local/extract_tiara_hits/main'
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

workflow ASCC_GENOMIC_ANALYSIS {
    take:
    reference_tuple_from_GG      // Channel: reference tuple from GENERATE_GENOME
    organellar_genomes           // Channel: tuple(meta, reference)
    fcs_db                       // [path(path)]
    reads                        // Channel: reads
    pacbio_database              // tuple [[meta.id], pacbio_database]
    ncbi_taxonomy_path
    ncbi_ranked_lineage_path
    nt_database_path
    diamond_nr_db_path
    diamond_uniprot_db_path
    taxid                        // NEW PARAMETER from dev branch
    nt_kraken_db_path           // NEW PARAMETER from dev branch

    main:
    ch_versions = Channel.empty()

    // Initialize output channels
    ch_kmers = Channel.of([[],[]])
    ch_tiara = Channel.of([[],[]])
    ch_nt_blast = Channel.of([[],[]])
    ch_blast_lineage = Channel.of([[],[]])
    ch_btk_format = Channel.of([[],[]])
    nr_full = Channel.of([[],[]])
    nr_hits = Channel.of([[],[]])
    un_full = Channel.of([[],[]])
    un_hits = Channel.of([[],[]])
    ch_mito = Channel.of([[],[]])
    ch_chloro = Channel.of([[],[]])
    ch_fcsadapt = Channel.of([[],[]])
    ch_fcsgx = Channel.of([[],[]])
    ch_coverage = Channel.of([[],[]])
    ch_bam = Channel.of([[],[]])
    ch_vecscreen = Channel.of([[],[]])
    ch_kraken1 = Channel.of([[],[]])
    ch_kraken2 = Channel.of([[],[]])
    ch_kraken3 = Channel.of([[],[]])
    ch_kmers_results = Channel.of([[],[]])

    //
    // KMERS ANALYSIS
    //
    if (params.run_kmers == "both" ||
        (params.run_kmers == "genomic" && params.genomic_only) ||
        (params.run_kmers == "organellar" && !params.genomic_only)) {
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

        ch_kmers_results    = GET_KMERS_PROFILE.out.kmers_results
                                .map { meta, dirs ->
                                    [[id: meta.id, process: "KMERS_RESULTS"], dirs]
                                }
                                .ifEmpty { [[],[]] }
    }

    //
    // TIARA ANALYSIS
    //
    if (params.run_tiara == "both" ||
        (params.run_tiara == "genomic" && params.genomic_only) ||
        (params.run_tiara == "organellar" && !params.genomic_only)) {
        EXTRACT_TIARA_HITS (
            reference_tuple_from_GG
        )
        ch_versions         = ch_versions.mix(EXTRACT_TIARA_HITS.out.versions)
        ch_tiara            = EXTRACT_TIARA_HITS.out.ch_tiara
                                .map { it ->
                                    [[id: it[0].id, process: "TIARA"], it[1]]
                                }
                                .ifEmpty { [[],[]] }
    }

    //
    // NT-BLAST ANALYSIS
    //
    if (params.run_nt_blast == "both" ||
        (params.run_nt_blast == "genomic" && params.genomic_only) ||
        (params.run_nt_blast == "organellar" && !params.genomic_only)) {
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
    }

    //
    // DIAMOND BLAST ANALYSIS
    //
    if (params.run_nr_diamond == "both" ||
        (params.run_nr_diamond == "genomic" && params.genomic_only) ||
        (params.run_nr_diamond == "organellar" && !params.genomic_only)) {
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
    }

    //
    // UNIPROT DIAMOND BLAST ANALYSIS
    //
    if (params.run_uniprot_diamond == "both" ||
        (params.run_uniprot_diamond == "genomic" && params.genomic_only) ||
        (params.run_uniprot_diamond == "organellar" && !params.genomic_only)) {
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
    }

    //
    // ORGANELLAR BLAST ANALYSIS
    //
    if (params.run_organellar_blast == "both" ||
        (params.run_organellar_blast == "genomic" && params.genomic_only) ||
        (params.run_organellar_blast == "organellar" && !params.genomic_only)) {
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
    }

    //
    // PACBIO BARCODE ANALYSIS
    //
    if (params.run_pacbio_barcodes == "both" ||
        (params.run_pacbio_barcodes == "genomic" && params.genomic_only) ||
        (params.run_pacbio_barcodes == "organellar" && !params.genomic_only)) {
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
    // FCS-ADAPTOR ANALYSIS
    //
    if (params.run_fcs_adaptor == "both" ||
        (params.run_fcs_adaptor == "genomic" && params.genomic_only) ||
        (params.run_fcs_adaptor == "organellar" && !params.genomic_only)) {
        RUN_FCSADAPTOR (
            reference_tuple_from_GG
        )
        ch_versions         = ch_versions.mix(RUN_FCSADAPTOR.out.versions)

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
    }

    //
    // FCS-GX ANALYSIS
    //
    if (params.run_fcsgx == "both" ||
        (params.run_fcsgx == "genomic" && params.genomic_only) ||
        (params.run_fcsgx == "organellar" && !params.genomic_only)) {
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

        ch_fcsgx            = RUN_FCSGX.out.fcsgxresult
                                .map { it ->
                                    [[id: it[0].id, process: "FCSGX result"], it[1]]
                                }
                                .ifEmpty { [[],[]] }
    }

    //
    // READ COVERAGE ANALYSIS
    //
    if (params.run_coverage == "both" ||
        (params.run_coverage == "genomic" && params.genomic_only) ||
        (params.run_coverage == "organellar" && !params.genomic_only)) {
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
    }

    //
    // VECTOR SCREENING ANALYSIS
    //
    if (params.run_vecscreen == "both" ||
        (params.run_vecscreen == "genomic" && params.genomic_only) ||
        (params.run_vecscreen == "organellar" && !params.genomic_only)) {
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
    }

    //
    // KRAKEN CLASSIFICATION
    //
    if (params.run_kraken == "both" ||
        (params.run_kraken == "genomic" && params.genomic_only) ||
        (params.run_kraken == "organellar" && !params.genomic_only)) {
        RUN_NT_KRAKEN(
            reference_tuple_from_GG,
            nt_kraken_db_path.first(),
            ncbi_ranked_lineage_path.first()
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
    }

    emit:
    versions = ch_versions
    kmers = ch_kmers
    tiara = ch_tiara
    nt_blast = ch_nt_blast
    blast_lineage = ch_blast_lineage
    btk_format = ch_btk_format
    nr_full = nr_full
    nr_hits = nr_hits
    un_full = un_full
    un_hits = un_hits
    fcsgx = ch_fcsgx
    coverage = ch_coverage
    bam = ch_bam
    vecscreen = ch_vecscreen
    kraken1 = ch_kraken1
    kraken2 = ch_kraken2
    kraken3 = ch_kraken3
    fcsadapt = ch_fcsadapt
    kmers_results = ch_kmers_results
    // Additional outputs needed by reporting workflow
    mito = ch_mito
    chloro = ch_chloro
    fcsgx_report_txt = (params.run_fcsgx == "both" ||
                        (params.run_fcsgx == "genomic" && params.genomic_only) ||
                        (params.run_fcsgx == "organellar" && !params.genomic_only)) ?
                        RUN_FCSGX.out.fcsgx_report : Channel.of([[],[]])
    fcsgx_taxonomy_rpt = (params.run_fcsgx == "both" ||
                        (params.run_fcsgx == "genomic" && params.genomic_only) ||
                        (params.run_fcsgx == "organellar" && !params.genomic_only)) ?
                        RUN_FCSGX.out.taxonomy_report : Channel.of([[],[]])
    barcode_check_filtered = (params.run_pacbio_barcodes == "both" ||
                            (params.run_pacbio_barcodes == "genomic" && params.genomic_only) ||
                            (params.run_pacbio_barcodes == "organellar" && !params.genomic_only)) ?
                        PACBIO_BARCODE_CHECK.out.filtered : Channel.of([[],[]])
}
