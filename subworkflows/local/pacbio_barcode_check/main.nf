//
// PACBIO_BARCODE_CHECK IDENTIFIED LOCATIONS OF BARCODE SEQUENCES IN THE INPUT ASSEMBLY
//

//
// MODULE IMPORT BLOCK
//
include { CHECK_BARCODE       } from '../../../modules/local/check/barcode/main'
include { BLAST_MAKEBLASTDB   } from '../../../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN        } from '../../../modules/nf-core/blast/blastn/main'
include { FILTER_BARCODE      } from '../../../modules/local/filter/barcode/main'

workflow PACBIO_BARCODE_CHECK {
    take:
    reference_tuple         // tuple    [[meta.id], reference ]
    barcode_names           // val      (csv-list-string)
    barcode_database        // tupe     [[meta.id], barcode_database]

    main:
    ch_versions     = Channel.empty()

    //
    // MODULE: RUN BLAST WITH GENOME AGAINST BARCODE DB
    //
    BLAST_BLASTN (
        reference_tuple,
        barcode_database
    )
    ch_versions     = ch_versions.mix(BLAST_BLASTN.out.versions)


    //
    // LOGIC: FOR I (MAPPED TO OTHER CHANNELS) IN CSV LIST RUN FILTER BLAST
    //
    Channel.of(barcode_names)
        .set {barcodes_name_list}

    barcodes_name_list
        .map { it ->
            tuple( it.split(',') )
        }
        .flatten()
        .set {barcode_list}


    //
    // MODULE: CREATE A FILTERED BLAST OUTPUT PER BARCODE
    //
    reference_tuple
        .combine(BLAST_BLASTN.out.txt, by:0)
        .combine(barcode_list)
        .set { filter_input }


    FILTER_BARCODE (
        filter_input
    )
    ch_versions     = ch_versions.mix(FILTER_BARCODE.out.versions)

    emit:
    filtered        = FILTER_BARCODE.out.debarcoded
    versions        = ch_versions
}
