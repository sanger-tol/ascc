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

workflow ASCC_ORGANELLAR {

    take:
    ch_samplesheet          // channel: samplesheet read in from --input
    validate_taxid_versions // Versions channel from main.nf
    include_steps           // params.include_steps
    exclude_steps           // params.exclude_steps
    fcs_db                  // path(file)
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
        "kmers", "tiara", "coverage", "nt_blast", "nr_diamond", "uniprot_diamond",
        "kraken", "fcs-gx", "fcs-adaptor", "vecscreen", "btk_busco",
        "pacbio_barcodes", "organellar_blast", "autofilter_assembly", "create_btk_dataset", "ALL", "NONE"]

    if (!full_list.containsAll(include_workflow_steps) && !full_list.containsAll(exclude_workflow_steps)) {
        exit 1, "There is an extra argument given on Command Line: \n Check contents of: $include_workflow_steps\nAnd $exclude_workflow_steps\nMaster list is: $full_list"
    }

    println "ORGANELLAR RUN -- INCLUDE STEPS INC.: $include_workflow_steps"
    println "ORGANELLAR RUN -- EXCLUDE STEPS INC.: $exclude_workflow_steps"


    //
    // LOGIC: CREATE btk_busco_run_mode VALUE
    //
    btk_busco_run_mode = params.btk_busco_run_mode ?: "conditional"


    //
    // LOGIC: PRETTY NOTIFICATION OF FILES AT STAGE
    //
    ch_samplesheet
        .map { meta, sample ->
            println "ORGANELLAR WORKFLOW:\n\t-- $meta\n\t-- $sample"
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
        ch_tiara            = EXTRACT_TIARA_HITS.out.ch_tiara.map{it[1]}
    } else {
        ch_tiara            = []
    }


    //
    // SUBWORKFLOW: IDENTITY PACBIO BARCODES IN INPUT DATA
    //
    if ( (include_workflow_steps.contains('pacbio_barcodes') || include_workflow_steps.contains('ALL')) && !exclude_workflow_steps.contains("pacbio_barcodes") ) {
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
    if ( (include_workflow_steps.contains('fcs-adaptor') || include_workflow_steps.contains('ALL')) && !exclude_workflow_steps.contains("fcs-adaptor") ) {
        RUN_FCSADAPTOR (
            ESSENTIAL_JOBS.out.reference_tuple_from_GG // Again should this be the validated fasta?
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
    if ( (include_workflow_steps.contains('fcs-gx') || include_workflow_steps.contains('ALL')) && !exclude_workflow_steps.contains("fcs-gx") ) {

        ESSENTIAL_JOBS.out.reference_tuple_from_GG
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

        ch_fcsgx            = RUN_FCSGX.out.fcsgxresult.map{it[1]}
        ch_versions         = ch_versions.mix(RUN_FCSGX.out.versions)
    } else {
        ch_fcsgx            = []
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
    // SUBWORKFLOW: SCREENING FOR VECTOR SEQUENCE
    //
    if ( (include_workflow_steps.contains('vecscreen') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("vecscreen")
    ) {
        RUN_VECSCREEN (
            ESSENTIAL_JOBS.out.reference_tuple_from_GG, // Again should this be the validated fasta?
            params.vecscreen_database_path
        )
        ch_vecscreen        = RUN_VECSCREEN.out.vecscreen_contam.map{it[1]}
        ch_versions         = ch_versions.mix(RUN_VECSCREEN.out.versions)
    } else {
        ch_vecscreen        = []
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
        ch_kraken1          = RUN_NT_KRAKEN.out.classified.map{it[1]}
        ch_kraken2          = RUN_NT_KRAKEN.out.report.map{it[1]}
        ch_kraken3          = RUN_NT_KRAKEN.out.lineage

        ch_versions         = ch_versions.mix(RUN_NT_KRAKEN.out.versions)
    } else {
        ch_kraken1          = []
        ch_kraken2          = []
        ch_kraken3          = []
    }


    //
    // LOGIC: WE NEED TO MAKE SURE THAT THE INPUT SEQUENCE IS OF AT LEAST LENGTH OF params.seqkit_window
    //
    ESSENTIAL_JOBS.out.reference_with_seqkit
        //
        // Here we are using the un-filtered genome, any filtering may (accidently) cause an empty fasta
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
        .set{ valid_length_fasta }

    valid_length_fasta
        .map{ meta, file ->
            println "Running BLAST (NT, DIAMOND, NR) on VALID ORGANELLE: $meta --- $file"
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
        // NOTE: ch_nt_blast needs to be set in two places incase it
        //          fails during the run (This IS an expected outcome of this subworkflow)
        //
        ch_nt_blast         = []
        ch_blast_lineage    = []


        SUBWORKFLOW: EXTRACT RESULTS HITS FROM NT-BLAST

        EXTRACT_NT_BLAST (
            valid_length_fasta,
            Channel.value(params.nt_database_path),
            Channel.value(params.ncbi_accession_ids_folder),
            Channel.value(params.ncbi_ranked_lineage_path)
        )
        ch_versions         = ch_versions.mix(EXTRACT_NT_BLAST.out.versions)
        ch_nt_blast         = EXTRACT_NT_BLAST.out.ch_blast_hits.map{it[1]}
        ch_blast_lineage    = EXTRACT_NT_BLAST.out.ch_top_lineages.map{it[1]}

    } else {
        ch_nt_blast         = Channel.empty()
        ch_blast_lineage    = Channel.empty()
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
        nr_full             = NR_DIAMOND.out.reformed.map{it[1]}
        nr_hits             = NR_DIAMOND.out.hits_file.map{it[1]}
        ch_versions         = ch_versions.mix(NR_DIAMOND.out.versions)
    } else {
        nr_hits             = []
        nr_full             = []
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
        un_full             = UP_DIAMOND.out.reformed.map{it[1]}
        un_hits             = UP_DIAMOND.out.hits_file.map{it[1]}
        ch_versions         = ch_versions.mix(UP_DIAMOND.out.versions)
    } else {
        un_hits             = []
        un_full             = []
    }


    if ( (include_workflow_steps.contains('create_btk_dataset') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("create_btk_dataset")
    ) {
        ch_dot_genome           = ESSENTIAL_JOBS.out.dot_genome.map{it[1]}

        CREATE_BTK_DATASET (
            ESSENTIAL_JOBS.out.reference_tuple_from_GG,
            ch_dot_genome,
            [], //ch_kmers
            ch_tiara,
            ch_nt_blast,
            [], //ch_fcsgx,
            ch_bam,
            ch_coverage,
            ch_kraken1,
            ch_kraken2,
            ch_kraken3,
            nr_full,
            un_full,
            Channel.fromPath(params.ncbi_taxonomy_path).first()
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
