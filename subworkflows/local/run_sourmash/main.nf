#!/usr/bin/env nextflow
//
// MODULE IMPORT BLOCK
//

include { SOURMASH_SKETCH       }   from '../../../modules/nf-core/sourmash/sketch/main'

// DEPENDING ON WHETHER YOU DEV IN ASCC OR GO FOR NF-CORE ACCEPTANCE FIRST
// include { SOURMASH_MULTISEARCH  }   from '../../../modules/nf-core/sourmash/multisearch/main'
include { SOURMASH_MULTISEARCH  }   from '../../../modules/local/sourmash/multisearch/main'

include { CAT_CAT               }   from '../../../modules/nf-core/cat/cat/main'

// DEPENDING ON WHETHER YOU DEV IN ASCC OR GO FOR SANGER_TOL ACCEPTANCE FIRST
include { PARSE_SOURMASH } from '../../../modules/local/parse/sourmash/main'
// include { PARSE_SOURMASH } from '../../../modules/sanger-tol/parse/sourmash/main'



workflow RUN_SOURMASH {
    take:
    genome_fasta         // val(meta), path(fasta)
    sourmash_database       // val([db_paths])
    assembly_taxa_db        // path(assembly_taxa_db)
    target_taxa             // val(target_taxa)
    k                       // val(ksize)
    s                       // val(scaled)


    main:
    ch_versions     = Channel.empty()

    //
    // MODULE: CREATE SOURMASH DB FOR INPUT ASSEMBLY
    //
    SOURMASH_SKETCH (
        genome_fasta
    )

    ch_versions     = ch_versions.mix( SOURMASH_SKETCH.out.versions )


    //
    // LOGIC: CREATE BALANCED CHANNELS
    //          (signatures WILL BE MULTIPLIED BY DB COUNT)
    //

    multisearch_input = SOURMASH_SKETCH.out.signatures.combine(sourmash_database) // [[meta], signature, database] for db in dbs

    //
    // MODULE: SOURMASH MULTISEARCH RUNS THE GENOME OVER X NUMBER OF INPUT
    //          DATABASES

    SOURMASH_MULTISEARCH (
        multisearch_input,
        k,
        s
    )
    ch_versions     = ch_versions.mix( SOURMASH_MULTISEARCH.out.versions )


    //
    // MODULE: MERGE MULTISEATCH RESULTS
    //

    cat_cat_input = SOURMASH_MULTISEARCH.out.multisearch_results.groupTuple().view()

    CAT_CAT (
        cat_cat_input
    )
    ch_versions     = ch_versions.mix( CAT_CAT.out.versions )


    //
    // MODULE: PARSE SOURMASH OUTPUT INTO USABLE FORM
    //

    PARSE_SOURMASH (
        CAT_CAT.out.file_out.view(),
        assembly_taxa_db,
        target_taxa
    )
    ch_versions     = ch_versions.mix( PARSE_SOURMASH.out.versions )

    emit:

    // sourmash_summary            = PARSE_SOURMASH.out.summary_csv
    // sourmash_non_target         = PARSE_SOURMASH.out.nontarget_csv
    versions                    = ch_versions
}
