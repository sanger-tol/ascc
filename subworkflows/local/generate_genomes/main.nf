//
// MODULE IMPORT BLOCK
//
include { SAMTOOLS_FAIDX            } from '../../../modules/nf-core/samtools/faidx/main'
include { GNU_SORT                  } from '../../../modules/nf-core/gnu/sort'

workflow GENERATE_GENOME {
    take:
    to_chromsize    // tuple [[meta.id], file]
    barcodes

    main:
    ch_versions     = channel.empty()

    //
    // MODULE: GENERATE INDEX OF REFERENCE
    //          EMITS REFERENCE INDEX FILE MODIFIED FOR SCAFF SIZES
    //
    SAMTOOLS_FAIDX (
        to_chromsize,
        channel.of([[],[]]),
        true
    )


    //
    // MODULE: SORT CHROM SIZES BY CHOM SIZE NOT NAME
    //
    GNU_SORT (
        SAMTOOLS_FAIDX.out.sizes
    )
    ch_versions     = ch_versions.mix(  GNU_SORT.out.versions )


    emit:
    dot_genome      = GNU_SORT.out.sorted
    ref_index       = SAMTOOLS_FAIDX.out.fai
    reference_tuple = to_chromsize
    versions        = ch_versions
}
