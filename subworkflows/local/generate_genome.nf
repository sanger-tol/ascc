#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//

workflow GENERATE_GENOME {
    take:
    taxid     // Channel val(taxid)
    reference  // Channel [ val(meta), path(file) ]

    main:
    ch_versions     = Channel.empty()

    //
    // LOGIC: GENERATES A REFERENCE DATA TUPLE
    //
    reference
        .combine( taxid )
        .map { it ->
            tuple ([taxid: it[1]],
                    it[0])
        }
        .set { reference_ch }

    emit:
    reference_tuple = reference_ch
    versions        = ch_versions.ifEmpty(null)
}