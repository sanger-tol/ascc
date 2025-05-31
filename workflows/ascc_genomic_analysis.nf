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
    ncbi_taxonomy_path
    ncbi_ranked_lineage_path
    nt_database_path
    diamond_nr_db_path
    diamond_uniprot_db_path
    taxid                        // NEW PARAMETER from dev branch
    nt_kraken_db_path           // NEW PARAMETER from dev branch

    main:
    ch_versions = Channel.empty()

    // Initialize empty channels for all outputs
    kmers = Channel.of([[],[]])
    tiara = Channel.of([[],[]])
    nt_blast = Channel.of([[],[]])
    blast_lineage = Channel.of([[],[]])
    btk_format = Channel.of([[],[]])
    nr_full = Channel.of([[],[]])
    nr_hits = Channel.of([[],[]])
    un_full = Channel.of([[],[]])
    un_hits = Channel.of([[],[]])
    fcsgx = Channel.of([[],[]])
    coverage = Channel.of([[],[]])
    bam = Channel.of([[],[]])
    vecscreen = Channel.of([[],[]])
    kraken1 = Channel.of([[],[]])
    kraken2 = Channel.of([[],[]])
    kraken3 = Channel.of([[],[]])
    fcsadapt = Channel.of([[],[]])
    kmers_results = Channel.of([[],[]])

    // Analysis workflows will be extracted from the monolithic file
    // This is a placeholder structure - content will be added in next step

    emit:
    versions = ch_versions
    kmers = kmers
    tiara = tiara
    nt_blast = nt_blast
    blast_lineage = blast_lineage
    btk_format = btk_format
    nr_full = nr_full
    nr_hits = nr_hits
    un_full = un_full
    un_hits = un_hits
    fcsgx = fcsgx
    coverage = coverage
    bam = bam
    vecscreen = vecscreen
    kraken1 = kraken1
    kraken2 = kraken2
    kraken3 = kraken3
    fcsadapt = fcsadapt
    kmers_results = kmers_results
}
