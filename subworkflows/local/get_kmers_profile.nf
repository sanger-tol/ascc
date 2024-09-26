#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { GET_KMER_COUNTS                       } from '../../modules/local/get_kmer_counts'
include { KMER_COUNT_DIM_REDUCTION              } from '../../modules/local/kmer_count_dim_reduction'
include { KMER_COUNT_DIM_REDUCTION_COMBINE_CSV  } from '../../modules/local/kmer_count_dim_reduction_combine_csv'

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
        // Needs to be combined here otherwise only a single channel will move forwards.
        .combine(kmer_size)
        .multiMap{ meta, file ,kmer ->
            assembly: tuple(meta, file)
            kmer_size: kmer
        }
        .set { modified_input }

    //
    // MODULE: PRODUCE KMER COUNTS (USING KCOUNTER)
    //
    GET_KMER_COUNTS (
        modified_input.assembly,      // val(meta), path(reads)
        modified_input.kmer_size      // val kmer_size
    )
    ch_versions = ch_versions.mix(GET_KMER_COUNTS.out.versions)

    //
    // LOGIC: CREATE CHANNEL OF LIST OF SELECTED METHODS
    //
    dimensionality_reduction_methods.splitCsv().flatten().set{ methods }
    autoencoder_epochs_count.flatten().unique().set{ epoch } // Why does this value end up as a queue channel? Check Yaml input!!!
    epoch.view{ " epochs: $it"}

    methods
        .combine(GET_KMER_COUNTS.out.csv)
        .combine(n_neighbors_setting)
        .combine(epoch)
        .multiMap { method, meta, files, neighbours, autoencoder ->
            csv: tuple([id: meta.id, single_end: true], files)
            methods: method
            n_neighbours: neighbours
            autoencoder_val: autoencoder
        }
        .set{dim_reduction}

    //
    // MODULE: DIMENSIONALITY REDUCTION OF KMER COUNTS, USING SPECIFIED METHODS
    //
    KMER_COUNT_DIM_REDUCTION (
        dim_reduction.csv,                // val(meta), path(kmer_counts)
        dim_reduction.methods,            // val dimensionality_reduction_method
        dim_reduction.n_neighbours,       // val n_neighbors_setting
        dim_reduction.autoencoder_val     // val autoencoder_epochs_count
    )
    ch_versions = ch_versions.mix(KMER_COUNT_DIM_REDUCTION.out.versions)

    //
    // LOGIC: PREPARING INPUT TO COMBINE OUTPUT CSV FOR EACH METHOD
    //
    KMER_COUNT_DIM_REDUCTION.out.csv
        .groupTuple()
        .set { collected_files_for_combine }

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
