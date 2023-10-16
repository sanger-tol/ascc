#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { GET_KMER_COUNTS } from '../../modules/local/get_kmers_count'
include { KMER_COUNT_DIM_REDUCTION } from '../../modules/local/kmer_count_dim_reduction'

workflow GET_KMERS_PROFILE {
    take:
    assembly_fasta                      // Channel [ val(meta), path(file) ]
    kmer_size                           // Channel [ val(integer) ]
    dimensionality_reduction_methods    // Channel [ val(string) ]
    n_neighbors_setting                 // Channel [ val(string) ]
    autoencoder_epochs_count            // Channel [ val(integer) ]

    main:
    ch_versions     = Channel.empty()

    assembly_fasta
        .map{ it ->
            tuple([id: it[0].id,
                    single_end: true
                ],
                it[1]
            )
        }
        .set { modified_input }

    //
    // MODULE: PRODUCE KMER COUNTS (USING KCOUNTER)
    //
    GET_KMER_COUNTS (
        modified_input,      // val(meta), path(reads)
        kmer_size            // val kmer_size
    )
    ch_versions = ch_versions.mix(GET_KMER_COUNTS.out.versions)

    //
    // MODULE: DIMENSIONALITY REDUCTION OF KMER COUNTS, USING SPECIFIED METHODS 
    //
    KMER_COUNT_DIM_REDUCTION (
        GET_KMER_COUNTS.out.csv,             // val(meta), path(kmer_counts)
        dimensionality_reduction_methods,    // val dimensionality_reduction_methods
        n_neighbors_setting,                 // val n_neighbors_setting
        autoencoder_epochs_count            // val autoencoder_epochs_count
    )
    ch_versions = ch_versions.mix(KMER_COUNT_DIM_REDUCTION.out.versions)

    emit:
    KMER_COUNT_DIM_REDUCTION.out.csv
    versions        = ch_versions.ifEmpty(null)
}
