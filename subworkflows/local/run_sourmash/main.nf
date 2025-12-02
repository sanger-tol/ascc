#!/usr/bin/env nextflow
//
// MODULE IMPORT BLOCK
//

// include { SOURMASH_SKETCH       }   from '../../../modules/nf-core/sourmash/sketch/main'
// Local version to pass arguments dynamically, they are defined depending on databases parameters (line 39-52)
include { SOURMASH_SKETCH       }   from '../../../modules/local/sourmash/sketch/main'

// DEPENDING ON WHETHER YOU DEV IN ASCC OR GO FOR NF-CORE ACCEPTANCE FIRST
// include { SOURMASH_MULTISEARCH  }   from '../../../modules/nf-core/sourmash/multisearch/main'
include { SOURMASH_MULTISEARCH  }   from '../../../modules/local/sourmash/multisearch/main'

include { CAT_CAT as CAT_MULTISEARCH_RESULTS }   from '../../../modules/nf-core/cat/cat/main'
include { CAT_CAT as CAT_TAXA_DB             }   from '../../../modules/nf-core/cat/cat/main'

// DEPENDING ON WHETHER YOU DEV IN ASCC OR GO FOR SANGER_TOL ACCEPTANCE FIRST
include { PARSE_SOURMASH        }   from '../../../modules/local/parse/sourmash/main'
// include { PARSE_SOURMASH        }   from '../../../modules/sanger-tol/parse/sourmash/main'

include { GET_TARGET_TAXA       }   from '../../../modules/local/get/target_taxa/main'



workflow RUN_SOURMASH {
    take:
    genome_fasta            // tuple val(meta), path(fasta) - meta should contain taxid
    sourmash_databases      // channel: list of database maps [name, path, k_available, k_for_search, s, assembly_taxa_db]
    ncbi_ranked_lineage     // path: NCBI ranked lineage file
    taxonomy_level          // val: taxonomy level for target_taxa extraction (e.g., 'order', 'family')


    main:
    ch_versions     = Channel.empty()

    //
    // LOGIC: COLLECT ALL DATABASE CONFIGURATIONS AND COMPUTE SKETCH PARAMETERS
    //
    // Collect all databases into a list to compute k and s parameters
    sourmash_databases
        .collect()
        .map { db_list ->
            // Compute sketch parameters from database configuration
            def uniqueK = db_list.collect { it.k_for_search }.unique().sort()
            def minS = db_list.collect { it.s }.min()
            def kParams = uniqueK.collect { "k=${it}" }.join(',')
            def sketch_args = "dna -p scaled=${minS},${kParams}"

            log.info "[RUN_SOURMASH] Computed sketch parameters: ${sketch_args}"
            [sketch_args, db_list]  // Return both sketch params and db_list
        }
        .set { ch_sketch_and_dbs }

    // Split into separate channels
    ch_sketch_params = ch_sketch_and_dbs.map { it[0] }
    ch_all_databases = ch_sketch_and_dbs.map { it[1] }

    //
    // MODULE: CREATE SOURMASH SKETCH FOR INPUT ASSEMBLY
    //
    // Sketch is created once with all k values needed across all databases
    // Pass dynamic sketch parameters computed from database configurations
    SOURMASH_SKETCH (
        genome_fasta,
        ch_sketch_params
    )
    ch_versions = ch_versions.mix(SOURMASH_SKETCH.out.versions)

    //
    // LOGIC: CREATE BALANCED CHANNELS FOR MULTISEARCH
    //
    // Combine signature with each database, pass k and s as separate parameters
    SOURMASH_SKETCH.out.signatures
        .combine(sourmash_databases)
        .map { meta, signature, db_config ->
            [meta, signature, file(db_config.path), db_config.k_for_search, db_config.s]
        }
        .set { ch_multisearch_input }

    //
    // MODULE: SOURMASH MULTISEARCH
    //
    // Run multisearch for each database with its specific k_for_search and s
    SOURMASH_MULTISEARCH (
        ch_multisearch_input
    )
    ch_versions = ch_versions.mix(SOURMASH_MULTISEARCH.out.versions)

    //
    // MODULE: MERGE MULTISEARCH RESULTS
    //
    // Group all results by sample id and concatenate
    SOURMASH_MULTISEARCH.out.multisearch_results
        .groupTuple()
        .set { ch_cat_multisearch_input }

    CAT_MULTISEARCH_RESULTS (
        ch_cat_multisearch_input
    )
    ch_versions = ch_versions.mix(CAT_MULTISEARCH_RESULTS.out.versions)

    //
    // LOGIC: PREPARE ASSEMBLY_TAXA_DB
    //
    // Collect all unique assembly_taxa_db paths and merge if needed
    ch_all_databases
        .map { db_list ->
            def taxa_paths = db_list.collect { it.assembly_taxa_db }.unique()

            // Check for duplicates
            def all_taxa_paths = db_list.collect { it.assembly_taxa_db }
            def duplicates = all_taxa_paths.findAll { p -> all_taxa_paths.count(p) > 1 }.unique()
            if (duplicates.size() > 0) {
                log.warn "[RUN_SOURMASH] Multiple databases share the same assembly_taxa_db: ${duplicates.join(', ')}"
            }

            taxa_paths
        }
        .flatten()
        .map { path -> file(path) }
        .collect()
        .set { ch_taxa_files }

    // Merge taxa files if multiple
    ch_taxa_files
        .map { files ->
            if (files.size() > 1) {
                [[id: 'merged_taxa'], files]
            } else {
                [[id: 'single_taxa'], files[0]]
            }
        }
        .set { ch_taxa_to_merge }

    CAT_TAXA_DB (
        ch_taxa_to_merge
    )
    ch_versions = ch_versions.mix(CAT_TAXA_DB.out.versions)

    //
    // MODULE: GET TARGET TAXA FROM TAXID
    //

    // Extract taxid from meta and get target_taxa
    genome_fasta
        .map { meta, fasta ->
            def taxid = meta.taxid ?: "UNKNOWN"
            [meta, taxid]
        }
        .set { ch_taxid_input }

    GET_TARGET_TAXA (
        ch_taxid_input,
        ncbi_ranked_lineage,
        taxonomy_level
    )
    ch_versions = ch_versions.mix(GET_TARGET_TAXA.out.versions)

    // Read target_taxa from file content and convert to value
    GET_TARGET_TAXA.out.target_taxa
        .map { meta, target_taxa_file ->
            def target_taxa_str = target_taxa_file.text.trim()
            [meta.id, target_taxa_str]
        }
        .set { ch_target_taxa_val }

    // Join merged multisearch results with target_taxa by sample id
    CAT_MULTISEARCH_RESULTS.out.file_out
        .map { meta, file -> [meta.id, meta, file] }
        .join(ch_target_taxa_val)
        .map { id, meta, results_file, target_taxa_str ->
            [meta, results_file, target_taxa_str]
        }
        .set { ch_results_with_target }

    //
    // MODULE: PARSE SOURMASH OUTPUT
    //
    // Extract assembly_taxa_db file and target_taxa value
    ch_taxa_db_file = CAT_TAXA_DB.out.file_out.map { meta, file -> file }
    ch_target_taxa_val = ch_results_with_target.map { meta, results, target_taxa -> target_taxa }

    PARSE_SOURMASH (
        ch_results_with_target.map { meta, results, target_taxa -> [meta, results] },
        ch_taxa_db_file,
        ch_target_taxa_val
    )
    ch_versions = ch_versions.mix(PARSE_SOURMASH.out.versions)

    emit:
    sourmash_summary            = PARSE_SOURMASH.out.multisearch_summary
    sourmash_non_target         = PARSE_SOURMASH.out.multisearch_non_target
    versions                    = ch_versions
}
