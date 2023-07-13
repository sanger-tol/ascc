#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { KRAKEN2_KRAKEN2 } from '../../modules/nf-core/kraken2/kraken2/main'
include { GET_LINEAGE_FOR_KRAKEN } from '../../modules/local/get_lineage_for_kraken'

workflow RUN_NT_KRAKEN {
    take:
    assembly_fasta
    nt_kraken_db_path
    ncbi_rankedlineage_path

    main:
    ch_versions     = Channel.empty()

    //
    // MODULE: Kraken2 run on assembly fasta.
    //
    KRAKEN2_KRAKEN2 ( 
        assembly_fasta,      // val(meta), path(reads)
        nt_kraken_db_path,   // path db
        false,               // val save_output_fastqs
        true                 // val save_reads_assignment
    )
    ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions)

    //
    // MODULE: Get lineage for kraken output.
    //
    GET_LINEAGE_FOR_KRAKEN ( 
        KRAKEN2_KRAKEN2.out.classified_reads_assignment,
        ncbi_rankedlineage_path
    )
    ch_versions = ch_versions.mix(GET_LINEAGE_FOR_KRAKEN.out.versions)
    
    emit:
    KRAKEN2_KRAKEN2.out.classified_reads_assignment
    KRAKEN2_KRAKEN2.out.report
    GET_LINEAGE_FOR_KRAKEN.out.txt
    versions        = ch_versions.ifEmpty(null)
}