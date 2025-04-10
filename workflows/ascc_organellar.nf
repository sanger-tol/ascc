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

include { paramsSummaryMultiqc                          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                        } from '../subworkflows/local/utils_nfcore_ascc_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASCC_ORGANELLAR {

    take:
    ch_samplesheet          // channel: samplesheet read in from --input
    include_steps_OBSELETE  // params.include_steps
    exclude_steps_OBSELETE  // params.exclude_steps
    fcs_db                  // path(file)
    reads
    scientific_name         // val(name)
    pacbio_database         // tuple [[meta.id], pacbio_database]

    main:
    ch_versions = Channel.empty()

    //
    // LOGIC: CONTROL OF THE INCLUDE AND EXCLUDE FLAGS
    //      TODO: THESE SHOULD CREATE A SET OF INCLUDE - EXCLUDE
    //      TODO: YES THIS IS DUPLICATED FROM PIPELINE INIT,
    //              HOWEVER THAT CONVERTED THE VALUES INTO A CHANNEL WHICH ISN'T THE EASIEST THING TO THEN PARSE OUT
    include_workflow_steps  = params.include ? params.include.split(",") : "ALL"
    exclude_workflow_steps  = params.exclude ? params.exclude.split(",") : "NONE"

    full_list               = [
        "kmers", "tiara", "coverage", "nt_blast", "nr_diamond", "uniprot_diamond",
        "kraken", "fcs-gx", "fcs-adaptor", "vecscreen", "btk_busco",
        "pacbio_barcodes", "organellar_blast", "autofilter_assembly", "create_btk_dataset", "ALL", "NONE"]

    if (!full_list.containsAll(include_workflow_steps) && !full_list.containsAll(exclude_workflow_steps)) {
        exit 1, "There is an extra argument given on Command Line: \n Check contents of: $include_workflow_steps\nAnd $exclude_workflow_steps\nMaster list is: $full_list"
    }

    log.info "ORGANELLAR RUN -- INCLUDE STEPS INC.: $include_workflow_steps"
    log.info "ORGANELLAR RUN -- EXCLUDE STEPS INC.: $exclude_workflow_steps"


    //
    // LOGIC: CREATE btk_busco_run_mode VALUE
    //
    btk_busco_run_mode = params.btk_busco_run_mode ?: "conditional"


    //
    // LOGIC: PRETTY NOTIFICATION OF FILES AT STAGE
    //
    ch_samplesheet
        .map { meta, sample ->
            log.info "ORGANELLAR WORKFLOW:\n\t-- $meta\n\t-- $sample"
        }


    //
    // SUBWORKFLOW: RUNS FILTER_FASTA, GENERATE .GENOME, CALCS GC_CONTENT AND FINDS RUNS OF N's
    //
    ESSENTIAL_JOBS(
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(ESSENTIAL_JOBS.out.versions)


    //
    // SUBWORKFLOW: EXTRACT RESULTS HITS FROM TIARA
    //
    if ( (include_workflow_steps.contains('tiara') || include_workflow_steps.contains('ALL')) && !exclude_workflow_steps.contains("tiara") ) {
        EXTRACT_TIARA_HITS (
            ESSENTIAL_JOBS.out.reference_tuple_from_GG
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
    if ( (include_workflow_steps.contains('pacbio_barcodes') || include_workflow_steps.contains('ALL')) && !exclude_workflow_steps.contains("pacbio_barcodes") ) {

        ESSENTIAL_JOBS.out.reference_tuple_from_GG
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
    if ( (include_workflow_steps.contains('fcs-adaptor') || include_workflow_steps.contains('ALL')) && !exclude_workflow_steps.contains("fcs-adaptor") ) {
        RUN_FCSADAPTOR (
            ESSENTIAL_JOBS.out.reference_tuple_from_GG
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
            // TODO: IS THIS AN ISSUE?

    } else {
        ch_fcsadapt         = Channel.of([[],[]])
    }


    //
    // SUBWORKFLOW: RUN FCS-GX TO IDENTIFY CONTAMINATION IN THE ASSEMBLY
    //
    if ( (include_workflow_steps.contains('fcs-gx') || include_workflow_steps.contains('ALL')) && !exclude_workflow_steps.contains("fcs-gx") ) {

        joint_channel = ESSENTIAL_JOBS.out.reference_tuple_from_GG
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
            ESSENTIAL_JOBS.out.reference_tuple_from_GG, // Again should this be the validated fasta?
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
            ESSENTIAL_JOBS.out.reference_tuple_from_GG, // Again should this be the validated fasta?
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
    // SUBWORKFLOW: RUN THE KRAKEN CLASSIFIER
    //
    if ( (include_workflow_steps.contains('kraken') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("kraken")
    ) {
        RUN_NT_KRAKEN(
            ESSENTIAL_JOBS.out.reference_tuple_from_GG,
            params.nt_kraken_database_path,
            params.ncbi_ranked_lineage_path
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
        ch_kraken1 = Channel.of([[],[]])
        ch_kraken2 = Channel.of([[],[]])
        ch_kraken3 = Channel.of([[],[]])

    }


    //
    // LOGIC: WE NEED TO MAKE SURE THAT THE INPUT SEQUENCE IS OF AT LEAST LENGTH OF params.seqkit_window
    //
    valid_length_fasta = ESSENTIAL_JOBS.out.reference_with_seqkit
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
            log.info "Running BLAST (NT, DIAMOND, NR) on VALID ORGANELLE: $meta --- $file"
        }

    //
    // LOGIC: THIS CONDITIONAL SHOULD EXECUTE THE PROCESS WHEN:
    //          INCLUDE STEPS ARE EITHER nt_blast AND all
    //              _AS WELL AS_
    //          EXCLUDE _NOT_ CONTAINING nt_blast AND THE valid_length_fasta IS NOT EMPTY
    //
    if ( (include_workflow_steps.contains('nt_blast') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("nt_blast") && !valid_length_fasta.ifEmpty(true)
    ) {

        //
        //SUBWORKFLOW: EXTRACT RESULTS HITS FROM NT-BLAST
        //

        EXTRACT_NT_BLAST (
            valid_length_fasta,
            Channel.value(params.nt_database_path),
            Channel.value(params.ncbi_accession_ids_folder),
            Channel.value(params.ncbi_ranked_lineage_path)
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

    } else {
        ch_nt_blast         = Channel.of( [[],[]] )
        ch_blast_lineage    = Channel.of( [[],[]] )
    }


    //
    // SUBWORKFLOW: DIAMOND BLAST FOR INPUT ASSEMBLY
    //
    if ( (include_workflow_steps.contains('nr_diamond') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("nr_diamond") && !valid_length_fasta.ifEmpty(true)
    ) {
        NR_DIAMOND (
            valid_length_fasta,
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
    //
    //qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms sphylums salltitles
    if ( (include_workflow_steps.contains('uniprot_diamond') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("uniprot_diamond") && !valid_length_fasta.ifEmpty(true)
    ) {
        UP_DIAMOND (
            valid_length_fasta,
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


    if ( (include_workflow_steps.contains('create_btk_dataset') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("create_btk_dataset")
    ) {

        //
        // LOGIC: FOUND RACE CONDITION EFFECTING LONG RUNNING JOBS
        //          AND INPUT TO HERE ARE NOW MERGED AND MAPPED
        //          EMPTY CHANNELS ARE CHECKED AND DEFAULTED TO [[],[]]
        //
        ch_organellar_cbtk_input = ESSENTIAL_JOBS.out.reference_tuple_from_GG
            .map{ it -> tuple([
                id: it[0].id,
                taxid: it[0].taxid,
                sci_name: it[0].sci_name,
                process: "REFERENCE"], it[1])
            }
            .mix(
                ESSENTIAL_JOBS.out.dot_genome.map{ it -> tuple([id: it[0].id, process: "GENOME"], it[1])},
                ch_tiara,
                ch_nt_blast,
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
            Channel.fromPath(params.ncbi_taxonomy_path).first(),
            scientific_name
        )
        ch_versions             = ch_versions.mix(CREATE_BTK_DATASET.out.versions)
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
