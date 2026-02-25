

include { SOURMASH_SKETCH                    }   from '../../../modules/local/sourmash/sketch/main'
include { SOURMASH_MULTISEARCH               }   from '../../../modules/local/sourmash/multisearch/main'
include { PARSE_SOURMASH                     }   from '../../../modules/local/parse/sourmash/main'
include { GET_TARGET_TAXA                    }   from '../../../modules/local/get/target_taxa/main'


workflow RUN_SOURMASH {
    take:
    reference_tuple         // channel: [meta, fasta] - meta.taxid required
    sourmash_databases      // channel: [name, path, k_available, k_for_search, s, assembly_taxa_db]
    ncbi_ranked_lineage     // path
    taxonomy_level          // val

    main:
    ch_versions     = Channel.empty()

    // Collect all databases and compute unified sketch parameters (k and scaled)
    sourmash_databases
        .collect()
        .map { db_list ->
            def uniqueK = db_list.collect { it.k_for_search }.unique().sort()
            def minS = db_list.collect { it.s }.min()
            def kParams = uniqueK.collect { "k=${it}" }.join(',')
            def sketch_args = "dna -p scaled=${minS},${kParams}"

            log.info "[RUN_SOURMASH] Computed sketch parameters: ${sketch_args}"
            [sketch_args, db_list]
        }
        .multiMap { sketch_args, db_list ->
            params: sketch_args
            databases: db_list
        }
        .set { ch_sketch_and_dbs }

    ch_sketch_params = ch_sketch_and_dbs.params
    ch_all_databases = ch_sketch_and_dbs.databases

    // Create sketch once with all k values needed across databases
    SOURMASH_SKETCH (
        reference_tuple,
        ch_sketch_params
    )
    ch_versions = ch_versions.mix(SOURMASH_SKETCH.out.versions)

    // Combine signature with each database (one search per database)
    SOURMASH_SKETCH.out.signatures
        .combine(sourmash_databases)
        .map { meta, signature, db_config ->
            [meta, signature, file(db_config.path), db_config.k_for_search, db_config.s]
        }
        .set { ch_multisearch_input }

    SOURMASH_MULTISEARCH (
        ch_multisearch_input
    )
    ch_versions = ch_versions.mix(SOURMASH_MULTISEARCH.out.versions)

    // Group multisearch results per sample
    SOURMASH_MULTISEARCH.out.multisearch_results
        .groupTuple()
        .set { ch_multisearch_results_grouped }

    // Collect unique assembly_taxa_db files from all databases
    ch_all_databases
        .map { db_list ->
            def taxa_paths = db_list.collect { it.assembly_taxa_db }.unique()
            taxa_paths
        }
        .flatten()
        .map { path -> file(path) }
        .collect()
        .set { ch_taxa_db_files }

    // Extract target taxa from NCBI taxonomy using meta.taxid
    GET_TARGET_TAXA (
        reference_tuple.map { meta, fasta -> [meta, meta.taxid ?: "UNKNOWN"] },
        ncbi_ranked_lineage,
        taxonomy_level
    )
    ch_versions = ch_versions.mix(GET_TARGET_TAXA.out.versions)


    ch_multisearch_results_grouped
        .map { meta, files -> [meta, files] }
        .join(GET_TARGET_TAXA.out.target_taxa.map { meta, tf -> [meta, tf.text.trim()] })
        .multiMap { meta, results_files, target_taxa_str ->
            results: tuple(meta, results_files)
            target_taxa: target_taxa_str
        }
        .set { ch_parse_input }

    // Debug: verify both channels have data
    ch_parse_input.results
        .view { meta, results -> "[DEBUG PARSE INPUT] id=${meta.id} num_results=${results.size()}" }

    ch_parse_input.target_taxa
        .view { target -> "[DEBUG TARGET TAXA] ${target}" }

    // Parse results and taxonomy database
    PARSE_SOURMASH (
        ch_parse_input.results,
        ch_taxa_db_files,
        ch_parse_input.target_taxa
    )
    ch_versions = ch_versions.mix(PARSE_SOURMASH.out.versions)

    emit:
    sourmash_summary            = PARSE_SOURMASH.out.multisearch_summary
    sourmash_non_target         = PARSE_SOURMASH.out.multisearch_non_target
    versions                    = ch_versions
}
