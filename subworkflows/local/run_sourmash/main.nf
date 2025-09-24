#!/usr/bin/env nextflow
//
// MODULE IMPORT BLOCK
//

include { SOURMASH_SKETCH       }   from '../../../modules/nf-core/sourmash/sketch/main'

// DEPENDING ON WHETHER YOU DEV IN ASCC OR GO FOR NF-CORE ACCEPTANCE FIRST
// include { SOURMASH_MULTISEARCH  }   from '../../../modules/nf-core/sourmash/multisearch/main'
// include { SOURMASH_MULTISEARCH  }   from '../../../modules/local/sourmash/multisearch/main'

include { CAT_CAT               }   from '../../../modules/nf-core/cat/cat/main'

// DEPENDING ON WHETHER YOU DEV IN ASCC OR GO FOR SANGER_TOL ACCEPTANCE FIRST
// include { PARSE_SOURMASH } from '../../../modules/local/parse/sourmash/main'
// include { PARSE_SOURMASH } from '../../../modules/sanger-tol/parse/sourmash/main'



workflow RUN_VECSCREEN {
    take:
    reference_tuple         // val(meta), path(fasta)
    sourmash_database       // val([db_paths])
    // target_taxa             // val(target_taxa)


    main:
    ch_versions     = Channel.empty()

    //
    // MODULE: CREATE SOURMASH DB FOR INPUT ASSEMBLY
    //
    SOURMASH_SKETCH (
        reference_tuple
    )
    ch_versions     = ch_versions.mix( SOURMASH_SKETCH.out.versions )


    //
    // LOGIC: CREATE BALANCED CHANNELS
    //          (signatures WILL BE MULTIPLIED BY DB COUNT)
    //
    // SOURMASH_SKETCH.out.signatures
    //  .combine(sourmash_database)
    //  .branch { sig, db ->
    //      signature: sig
    //      database: db
    //  }
    //  .set { multisearch_input }

    //
    // MODULE: SOURMASH MULTISEARCH RUNS THE GENOME OVER X NUMBER OF INPUT
    //          DATABASES
    // SOURMASH_MULTISEARCH (
    //     multisearch_input.signature,
    //     multisearch_input.database
    // )
    // ch_versions     = ch_versions.mix( SOURMASH_MULTISEARCH.out.versions )


    //
    // MODULE: MERGE MULTISEATCH RESULTS
    //
    // cat_cat_input = SOURMASH_MULTISEARCH.out.multisearch_results.collect()

    // CAT_CAT (
    //     cat_cat_input
    // )
    // ch_versions     = ch_versions.mix( CAT_CAT.out.versions )


    //
    // MODULE: PARSE SOURMASH OUTPUT INTO USABLE FORM
    //
    // PARSE_SOURMASH (
    //     CAT_CAT.out.file_out,
    //     sourmash_db_assemblies_with_taxa.csv,
    //     target_taxa
    // )
    // ch_versions     = ch_versions.mix( PARSE_SOURMASH.out.versions )

    emit:

    // sourmash_summary            = PARSE_SOURMASH.out.summary_csv
    // sourmash_non_target         = PARSE_SOURMASH.out.nontarget_csv
    versions                    = ch_versions
}
