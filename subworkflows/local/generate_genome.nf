#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//

workflow GENERATE_GENOME {
    take:
    assembly_title     // Channel val(assembly_title)
    reference  // Channel [ val(meta), path(file) ]

    main:
    ch_versions     = Channel.empty()

    //
    // LOGIC: GENERATES A REFERENCE DATA TUPLE
    //
    reference
        .combine( assembly_title )
        .map { it ->
            tuple ([id: it[1]],
                    it[0])
        }
        .set { reference_ch }

    // THIS IS HERE FOR FUTURE EXPANSION

    emit:
    reference_tuple = reference_ch
    versions        = ch_versions.ifEmpty(null)
}
