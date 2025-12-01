//
// MODULE IMPORT BLOCK
//
include { CHECK_BARCODE     } from '../../../modules/local/check/barcode/main'
include { BLAST_MAKEBLASTDB } from '../../../modules/nf-core/blast/makeblastdb/main'


workflow PREPARE_BLASTDB {
    take:
    sample_id               // val      params.sample-id
    pacbio_data             // val      params.pacbio_data
    pacbio_type             // val      params.pacbio_type
    barcodes_file           // val      params.barcodes_file
    barcode_names           // val      (csv-list-string)

    main:
    ch_versions             = Channel.empty()


    //
    // MODULE: CHECK FOR KNOWN BARCODES IN SAMPLE DATA
    //
    CHECK_BARCODE (
        [[id: sample_id], pacbio_data],
        barcodes_file,
        barcode_names
    )
    ch_versions             = ch_versions.mix(CHECK_BARCODE.out.versions)


    //
    // LOGIC: ENSURE THE VALID CHANNEL IS MIXED WITH THE BARCODES CHANNEL
    //          ACTS AS A GATEKEEPER FOR THE FLOW
    //
    CHECK_BARCODE.out.result
        .filter{it.contains("barcodes")} // Indicates it is a valid barcode
        .combine( barcodes_file )
        .map {str_info, file ->
            tuple(
                [id: "BARCODE_TO_MAKEDB", info: str_info],
                file
            )
        }
        .set {ch_new_barcodes}


    //
    // MODULE: GENERATE BLAST DB ON PACBIO BARCODES
    //
    BLAST_MAKEBLASTDB (
        ch_new_barcodes
    )
    ch_versions             = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)


    emit:
    barcodes_blast_db       = BLAST_MAKEBLASTDB.out.db
    versions                = ch_versions

}
