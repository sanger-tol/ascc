#!/usr/bin/env nextflow
// MODULE IMPORT BLOCK

include { CHUNK_ASSEMBLY_FOR_VECSCREEN  }   from '../../modules/local/chunk_assembly_for_vecscreen'
include { NCBITOOLS_VECSCREEN           }   from '../../modules/nf-core/ncbitools/vecscreen/main'
include { FILTER_VECSCREEN_RESULTS      }   from '../../modules/local/filter_vecscreen_results'
include { SUMMARISE_VECSCREEN_OUTPUT    }   from '../../modules/local/summarise_vecscreen_output'


workflow RUN_VECSCREEN {
    take:
    reference_tuple             // val(meta), path(fasta)
    vecscreen_database_tuple    // val(db_path)

    main:
    ch_versions                 = Channel.empty()


    //
    // MODULE: CHUNKS THE ASSEMBLY INTO PIECES WITH FIXED LENGTH
    //
    CHUNK_ASSEMBLY_FOR_VECSCREEN(
        reference_tuple
    )
    ch_versions                 = ch_versions.mix( CHUNK_ASSEMBLY_FOR_VECSCREEN.out.versions )


    //
    // MODULE: RUNS NCBI VECSCREEN
    //
    NCBITOOLS_VECSCREEN(
        CHUNK_ASSEMBLY_FOR_VECSCREEN.out.chunked_assembly,
        [[id: "db"], vecscreen_database_tuple]
    )
    ch_versions                 = ch_versions.mix( NCBITOOLS_VECSCREEN.out.versions )


    //
    // MODULE: REFORMATS VECSCREEN OUTPUT AND FILTERS IT TO REMOVE NO-HIT RESULTS AND KEEP
    //          ONLY HITS
    //
    FILTER_VECSCREEN_RESULTS(
        NCBITOOLS_VECSCREEN.out.vecscreen_output
    )
    ch_versions                 = ch_versions.mix( FILTER_VECSCREEN_RESULTS.out.versions )


    //
    // MODULE: CONVERTS COORDINATES IN ASSEMBLY CHUNKS BACK TO COORDINATES IN THE WHOLE
    //          ASSEMBLY AND WRITES A REPORT FILE
    //
    SUMMARISE_VECSCREEN_OUTPUT(
        FILTER_VECSCREEN_RESULTS.out.filtered_vecscreen_outfile
    )
    ch_versions                 = ch_versions.mix( SUMMARISE_VECSCREEN_OUTPUT.out.versions )

    emit:
    vecscreen_contam            = SUMMARISE_VECSCREEN_OUTPUT.out.vecscreen_contamination
    versions                    = ch_versions.ifEmpty( null ) // channel: [ versions.yml ]
}
