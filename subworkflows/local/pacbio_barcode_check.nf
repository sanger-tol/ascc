import { CHECK_BARCODE          } from '../modules/local/check_barcode'
import { BLAST_MAKEBLASTDB      } from '../modules/local/blast/makeblastdb'
import { BLAST_BLASTN           } from '../modules/local/blast/blastn'
import { FILTER_BARCODE         } from '../modules/local/filter_barcode'

workflow PACBIO_BARCODE_CHECK () {
    take:
    reference_tuple
    pacbio_tuple
    barcode_file
    barcode_multiplex

    main:
    ch_versions             = Channel.empty()

    //
    // MODULE:
    //
    CHECK_BARCODE (
        pacbio_tuple
        barcode_file,
        barcode_multiplex
    )
    ch_versions     = ch_versions.mix(CHECK_BARCODE.out.versions)

    //
    // MODULE: GENERATE BLAST DB ON ORGANELLAR GENOME
    //
    BLAST_MAKEBLASTDB (
        barcode_file
    )
    ch_versions     = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

    //
    // MODULE: RUN BLAST WITH GENOME AGAINST ORGANELLAR GENOME
    //
    BLAST_BLASTN (
        reference_tuple,
        BLAST_MAKEBLASTDB.out.db
    )
    ch_versions     = ch_versions.mix(BLAST_BLASTN.out.versions)

    //
    // LOGIC: FOR I IN CSV LIST RUN FILTER BLAST
    //
    // TODO: CLAFFIFY THIS BIT

    FILTER_BARCODE (
        reference_tuple,
        //i in csv
    )
    ch_versions     = ch_versions.mix(FILTER_BARCODE.out.versions)

    emit:
    filtered        = FILTER_BARCODE.out.debarcoded
}
