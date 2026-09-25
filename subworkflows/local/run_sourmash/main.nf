

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
    exclude_accessions_file // path or [] - accessions whose hits are ignored (assembly_accession column)
    contig_sizes            // channel: [meta, sizes] (name<TAB>length, e.g. ESSENTIAL_JOBS dot_genome); joined on meta.id
    taxdump                 // path or [] - NCBI new_taxdump dir (sourmash_div, own-species exclusion)
    div_mapping             // path or [] - BLAST-div -> FCS-GX div table

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

    // Split reference_tuple explicitly so SKETCH and GET_TARGET_TAXA
    // each get an independent channel subscription (avoids implicit value-channel
    // consumption issues when ch_sketch_params is a single-emission channel).
    reference_tuple
        .multiMap { meta, fasta ->
            for_sketch:      [meta, fasta]
            for_target_taxa: [meta, meta.taxid ?: "UNKNOWN"]
        }
        .set { ch_ref_split }

    // Create sketch once with all k values needed across databases
    SOURMASH_SKETCH (
        ch_ref_split.for_sketch,
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
    // ncbi_ranked_lineage is a queue channel (single path); .first() converts it to
    // a value channel so GET_TARGET_TAXA can run once per sample (not just the first).
    GET_TARGET_TAXA (
        ch_ref_split.for_target_taxa,
        ncbi_ranked_lineage.first(),
        taxonomy_level
    )
    ch_versions = ch_versions.mix(GET_TARGET_TAXA.out.versions)


    // contig sizes are keyed on meta.id only (dot_genome carries [id:] meta); missing -> []
    ch_multisearch_results_grouped
        .join(GET_TARGET_TAXA.out.target_taxa.map { meta, tf -> [meta, tf.text.trim()] })
        .map { meta, results, target_taxa -> [meta.id, meta, results, target_taxa] }
        .join(contig_sizes.map { meta, sizes -> [meta.id, sizes] }, remainder: true)
        .filter { it.size() == 5 && it[1] != null }   // drop sizes-only entries (samples without results)
        .map { it -> [it[1], it[2], it[3], it[4] ?: []] }
        .set { ch_parse_input }

    // Parse results and taxonomy database
    PARSE_SOURMASH (
        ch_parse_input,
        ch_taxa_db_files,
        exclude_accessions_file,
        taxdump,
        div_mapping
    )
    ch_versions = ch_versions.mix(PARSE_SOURMASH.out.versions)

    emit:
    sourmash_summary            = PARSE_SOURMASH.out.multisearch_summary
    sourmash_non_target         = PARSE_SOURMASH.out.multisearch_non_target
    sourmash_classification     = PARSE_SOURMASH.out.multisearch_classification
    versions                    = ch_versions
}
