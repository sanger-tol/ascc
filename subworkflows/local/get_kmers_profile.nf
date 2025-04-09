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
    // Extract only the kmers_dim_reduction_embeddings.csv files from each directory
    KMER_COUNT_DIM_REDUCTION.out.csv
        .filter{meta, file -> !file.toString().contains("EMPTY")}
        .map { meta, files ->
            log.debug "Processing files for ${meta.id}: ${files}"
            def csv_files = []
            files.each { file ->
                log.debug "Checking file: ${file}"
                if (file.toString().endsWith("kmers_dim_reduction_embeddings.csv")) {
                    log.debug "Adding CSV file: ${file}"
                    csv_files.add(file)
                }
            }
            log.debug "Found ${csv_files.size()} CSV files for ${meta.id}"
            if (csv_files.size() > 0) {
                return tuple(meta, csv_files)
            } else {
                log.warn "No CSV files found for ${meta.id}, skipping"
                return null
            }
        }
        .filter { it != null }
        .groupTuple(by: [0])
        .map { meta, files ->
            log.debug "After groupTuple: ${meta.id} has ${files.size()} file groups"
            log.debug "Files: ${files.flatten()}"
            return tuple(meta, files.flatten())
        }
        .set { collected_files_for_combine }

    // Debug output using debug level logging
    collected_files_for_combine.map { meta, files ->
        log.debug "Files for combine: ${meta.id} -> ${files}"
        return tuple(meta, files)
    }
    .set { collected_files_for_combine }

    //
    // MODULE: COMBINE OUTPUTS OF MULTIPLE METHODS
    //
    KMER_COUNT_DIM_REDUCTION_COMBINE_CSV (
        collected_files_for_combine
    )
    ch_versions = ch_versions.mix(KMER_COUNT_DIM_REDUCTION_COMBINE_CSV.out.versions)

    // Collect the results directories from KMER_COUNT_DIM_REDUCTION
    KMER_COUNT_DIM_REDUCTION.out.results_dir
        .filter{meta, dir -> !dir.toString().contains("EMPTY")}
        .groupTuple(by: [0])
        .set { collected_results_dirs }

    emit:
    combined_csv = KMER_COUNT_DIM_REDUCTION_COMBINE_CSV.out.csv
    kmers_results = collected_results_dirs
    versions     = ch_versions.ifEmpty(null)
}
