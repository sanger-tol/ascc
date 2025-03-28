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
    // LOGIC: MODIFY THE INPUT TUPLE TO INCLUDE THE single_end VALUE
    //
    assembly_fasta
        .map{ meta, file ->
            tuple([id: meta.id,
                    single_end: true
                ],
                file
            )
        }
        .set { modified_input }


    //
    // MODULE: RUN KRAKEN2 ON THE INPUT GENOME - GENERATES TAXONOMIC CLASSIFICATION.
    //
    KRAKEN2_KRAKEN2 (
        modified_input,      // val(meta), path(reads)
        nt_kraken_db_path,   // path db
        false,               // val save_output_fastqs
        true                 // val save_reads_assignment
    )
    ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions)


    //
    // MODULE: GET LINEAGE FOR THE KRAKEN OUTPUT.
    //
    GET_LINEAGE_FOR_KRAKEN (
        KRAKEN2_KRAKEN2.out.classified_reads_assignment,
        ncbi_rankedlineage_path
    )
    ch_versions = ch_versions.mix(GET_LINEAGE_FOR_KRAKEN.out.versions)

    emit:
    classified      = KRAKEN2_KRAKEN2.out.classified_reads_assignment
    report          = KRAKEN2_KRAKEN2.out.report
    lineage         = GET_LINEAGE_FOR_KRAKEN.out.txt
    versions        = ch_versions
}
