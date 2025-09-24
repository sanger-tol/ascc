#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { SAMTOOLS_FAIDX            } from '../../../modules/nf-core/samtools/faidx/main'
include { GNU_SORT                  } from '../../../modules/nf-core/gnu/sort'
include { GET_LARGEST_SCAFFOLD      } from '../../../modules/local/get/largest_scaffold/main'

workflow GENERATE_GENOME {
    take:
    to_chromsize    // tuple [[meta.id], file]
    barcodes

    main:
    ch_versions     = Channel.empty()

    //
    // MODULE: GENERATE INDEX OF REFERENCE
    //          EMITS REFERENCE INDEX FILE MODIFIED FOR SCAFF SIZES
    //
    SAMTOOLS_FAIDX (
        to_chromsize,
        Channel.of([[],[]]),
        true
    )
    ch_versions     = ch_versions.mix(  SAMTOOLS_FAIDX.out.versions )


    //
    // MODULE: SORT CHROM SIZES BY CHOM SIZE NOT NAME
    //
    GNU_SORT (
        SAMTOOLS_FAIDX.out.sizes
    )
    ch_versions     = ch_versions.mix(  GNU_SORT.out.versions )


    //
    // MODULE: CUT OUT THE LARGEST SCAFFOLD SIZE AND USE AS A COMPARATOR AGAINST 512MB
    //          THIS IS THE CUT OFF FOR TABIX USING TBI INDEXES
    //
    GET_LARGEST_SCAFFOLD (
        SAMTOOLS_FAIDX.out.sizes
    )
    ch_versions     = ch_versions.mix( GET_LARGEST_SCAFFOLD.out.versions )

    emit:
    max_scaff_size  = GET_LARGEST_SCAFFOLD.out.scaff_size.toInteger()
    dot_genome      = GNU_SORT.out.sorted
    ref_index       = SAMTOOLS_FAIDX.out.fai
    reference_tuple = to_chromsize
    versions        = ch_versions
}
