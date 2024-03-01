#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { GET_KMER_COUNTS } from '../../modules/local/get_kmer_counts'
include { KMER_COUNT_DIM_REDUCTION } from '../../modules/local/kmer_count_dim_reduction'
include { KMER_COUNT_DIM_REDUCTION_COMBINE_CSV } from '../../modules/local/kmer_count_dim_reduction_combine_csv'

workflow GET_KMERS_PROFILE {
    take:
    assembly_fasta                      // Channel [ val(meta), path(file) ]
    kmer_size                           // Channel [ val(integer) ]
    dimensionality_reduction_methods    // Channel [ val(string) ]
    n_neighbors_setting                 // Channel [ val(string) ]
    autoencoder_epochs_count            // Channel [ val(integer) ]

    main:
    ch_versions     = Channel.empty()

    //
    // LOGIC: REFACTORING REFERENCE TUPLE
    //
    assembly_fasta
        .map{ meta, file ->
            tuple([id: meta.id,
                    single_end: true],
                file
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
    // LOGIC: CREATE CHANNEL OF LIST OF SELECTED METHODS
    //
    dimensionality_reduction_methods
        .splitCsv(sep: ',')
        .flatten()
        .combine( GET_KMER_COUNTS.out.csv )
        .combine( n_neighbors_setting )
        .combine( autoencoder_epochs_count )
        .multiMap { dim_red_method, csv_meta, csv_val, nn_setting, ae_setting ->
            kmer_csv :   tuple ([ id: csv_meta.id, single_end: true], csv_val)
            dimensionality_reduction_method : dim_red_method
            n_neighbors_setting : nn_setting
            autoencoder_epochs_count : ae_setting
        }
        .set{ ch_methods }

    //
    // MODULE: DIMENSIONALITY REDUCTION OF KMER COUNTS, USING SPECIFIED METHODS
    //
    KMER_COUNT_DIM_REDUCTION (
        ch_methods.kmer_csv,                         // val(meta), path(kmer_counts)
        ch_methods.dimensionality_reduction_method,  // val dimensionality_reduction_method
        ch_methods.n_neighbors_setting,              // val n_neighbors_setting
        ch_methods.autoencoder_epochs_count          // val autoencoder_epochs_count
    )
    ch_versions = ch_versions.mix(KMER_COUNT_DIM_REDUCTION.out.versions)

    //
    // LOGIC: PREPARING INPUT TO COMBINE OUTPUT CSV FOR EACH METHOD
    //
    KMER_COUNT_DIM_REDUCTION.out.csv
        .collect()
        .map { file ->
            tuple (
                    [   id  :   file[0].toString().split('/')[-1].split('_')[0] ], // Change sample ID
                    file
            )
        }
        .set { collected_files_for_combine }

    collected_files_for_combine.view()

    //
    // MODULE: COMBINE OUTPUTS OF MULTIPLE METHODS
    //
    KMER_COUNT_DIM_REDUCTION_COMBINE_CSV (
        collected_files_for_combine
    )
    ch_versions = ch_versions.mix(KMER_COUNT_DIM_REDUCTION_COMBINE_CSV.out.versions)

    emit:
    combined_csv = KMER_COUNT_DIM_REDUCTION_COMBINE_CSV.out.csv
    versions     = ch_versions.ifEmpty(null)
}
