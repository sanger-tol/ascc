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
    barcodes                = Channel.empty()

    if (barcode_file.isEmpty("YES") == "YES") {
        Channel
            .fromPath("./assets/pacbio_adaptors.fa")
            .map { it ->
                tuple(  [id: "pacbio_barcodes"],
                        it
                )
            }
            .set { barcodes }
    } else {
        Channel
            .fromPath(barcode_file)
            .map { it ->
                tuple(  [id: "pacbio_barcodes"],
                        it
                )
            }
            .set { barcodes }
    }


    barcodes.view()

    //
    // MODULE: CHECK FOR KNOWN BARCODES IN SAMPLE DATA
    //
    CHECK_BARCODE (
        pacbio_tuple
        barcodes,
        barcode_multiplex
    )
    ch_versions     = ch_versions.mix(CHECK_BARCODE.out.versions)

    //
    // MODULE: GENERATE BLAST DB ON ORGANELLAR GENOME
    //
    BLAST_MAKEBLASTDB (
        barcodes
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
