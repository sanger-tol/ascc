//
// MODULE IMPORT BLOCK
//
include { GET_KMER_COUNTS                       } from '../../../modules/local/get/kmer_counts/main'
include { KMER_COUNT_DIM_REDUCTION              } from '../../../modules/local/kmer_count/dim_reduction/main'
include { KMER_COUNT_DIM_REDUCTION_COMBINE_CSV  } from '../../../modules/local/kmer_count/dim_reduction_combine_csv/main'

workflow GET_KMERS_PROFILE {
    take:
    assembly_fasta                      // Channel [ val(meta), path(file) ]
    kmer_size                           // Channel [ val(integer) ]
    dimensionality_reduction_methods    // Channel [ val(string) ]
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
    Channel.fromList(params.dimensionality_reduction_methods)
        .set{dim_methods}

    Channel.from(params.n_neighbours)
        .set{hey_neighbour}

    dim_methods
        .combine(GET_KMER_COUNTS.out.csv)
        .combine(autoencoder_epochs_count.first())
        .combine(hey_neighbour)
        .multiMap { dr_method, csv_meta, csv_file, epochs, n_neighbours ->
            method_name: dr_method
            kmer_csv: [csv_meta, csv_file]
            epoch_count: epochs
            n_neighbors_setting: n_neighbours

        }
        .set{ dim_reduction }

    //
    // MODULE: DIMENSIONALITY REDUCTION OF KMER COUNTS, USING SPECIFIED METHODS
    //
    KMER_COUNT_DIM_REDUCTION (
        dim_reduction.kmer_csv,             // val(meta), path(kmer_counts)
        dim_reduction.method_name,          // val dimensionality_reduction_method
        dim_reduction.n_neighbors_setting,  // val n_neighbors_setting
        dim_reduction.epoch_count           // val autoencoder_epochs_count
    )
    ch_versions = ch_versions.mix(KMER_COUNT_DIM_REDUCTION.out.versions)

    //
    // LOGIC: PREPARING INPUT TO COMBINE OUTPUT CSV FOR EACH METHOD
    //
    KMER_COUNT_DIM_REDUCTION.out.kmers_dim_reduction_dir
        .filter{meta, file -> !file.toString().contains("EMPTY")}
        .groupTuple(by: [0])
        .set { collected_files_for_combine }

    // Collect the results directories from KMER_COUNT_DIM_REDUCTION
    KMER_COUNT_DIM_REDUCTION.out.results_dir
        .filter{meta, dir -> !dir.toString().contains("EMPTY")}
        .groupTuple(by: [0])
        .set { collected_results_dirs }

    //
    // MODULE: COMBINE OUTPUTS OF MULTIPLE METHODS
    //
    KMER_COUNT_DIM_REDUCTION_COMBINE_CSV (
        collected_files_for_combine
    )
    ch_versions = ch_versions.mix(KMER_COUNT_DIM_REDUCTION_COMBINE_CSV.out.versions)

    emit:
    combined_csv = KMER_COUNT_DIM_REDUCTION_COMBINE_CSV.out.csv
    kmers_results = collected_files_for_combine
    versions     = ch_versions
}
