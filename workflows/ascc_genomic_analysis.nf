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
    include_workflow_steps       // List of included workflow steps
    exclude_workflow_steps       // List of excluded workflow steps
    organellar_genomes           // Channel: tuple(meta, reference)
    fcs_db                       // [path(path)]
    reads                        // Channel: reads
    pacbio_database              // tuple [[meta.id], pacbio_database]

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

    //
    // KMERS ANALYSIS
    //
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
    }

    //
    // NT-BLAST ANALYSIS
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
    if ( (include_workflow_steps.contains('nr_diamond') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("nr_diamond")
    ) {
        NR_DIAMOND (
            reference_tuple_from_GG,
            params.diamond_nr_database_path
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
    if ( (include_workflow_steps.contains('uniprot_diamond') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("uniprot_diamond")
    ) {
        UP_DIAMOND (
            reference_tuple_from_GG,
            params.diamond_uniprot_database_path
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
    if ( (include_workflow_steps.contains('organellar_blast') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("organellar_blast")
    ) {
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
    if ( (include_workflow_steps.contains('pacbio_barcodes') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("pacbio_barcodes")
    ) {
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
    if ( (include_workflow_steps.contains('fcs-adaptor') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("fcs-adaptor")
    ) {
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
    }

    //
    // FCS-GX ANALYSIS
    //
    if ( (include_workflow_steps.contains('fcs-gx') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("fcs-gx")
    ) {
        joint_channel = reference_tuple_from_GG
            .combine(fcs_db)
            .combine(Channel.of(params.taxid))
            .combine(Channel.of(params.ncbi_ranked_lineage_path))
            .multiMap { meta, ref, db, taxid, tax_path ->
                reference: [meta, taxid, ref]
                fcs_db_path: db
                taxid_val: taxid
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
    }

    //
    // READ COVERAGE ANALYSIS
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
    }

    //
    // VECTOR SCREENING ANALYSIS
    //
    if ( (include_workflow_steps.contains('vecscreen') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("vecscreen")
    ) {
        RUN_VECSCREEN (
            reference_tuple_from_GG,
            params.vecscreen_database_path
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
    }

    //
    // KRAKEN CLASSIFICATION
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
    }

    emit:
    kmers = ch_kmers
    kmers_results = ch_kmers_results
    tiara = ch_tiara
    nt_blast = ch_nt_blast
    blast_lineage = ch_blast_lineage
    btk_format = ch_btk_format
    nr_full = nr_full
    nr_hits = nr_hits
    un_full = un_full
    un_hits = un_hits
    mito = ch_mito
    chloro = ch_chloro
    fcsadapt = ch_fcsadapt
    fcsgx = ch_fcsgx
    fcsgx_report_txt = include_workflow_steps.contains('fcs-gx') || include_workflow_steps.contains('ALL') ? 
                        RUN_FCSGX.out.fcsgx_report : Channel.of([[],[]])
    fcsgx_taxonomy_rpt = include_workflow_steps.contains('fcs-gx') || include_workflow_steps.contains('ALL') ? 
                        RUN_FCSGX.out.taxonomy_report : Channel.of([[],[]])
    coverage = ch_coverage
    bam = ch_bam
    vecscreen = ch_vecscreen
    kraken1 = ch_kraken1
    kraken2 = ch_kraken2
    kraken3 = ch_kraken3
    versions = ch_versions
}
