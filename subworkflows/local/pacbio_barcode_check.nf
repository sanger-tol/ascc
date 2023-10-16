//
// PACBIO_BARCODE_CHECK IDENTIFIED LOCATIONS OF BARCODE SEQUENCES IN THE INPUT ASSEMBLY
//

//
// MODULE IMPORT BLOCK
//
include { CHECK_BARCODE          } from '../../modules/local/check_barcode'
include { BLAST_MAKEBLASTDB      } from '../../modules/nf-core/blast/makeblastdb'
include { BLAST_BLASTN           } from '../../modules/nf-core/blast/blastn'
include { FILTER_BARCODE         } from '../../modules/local/filter_barcode'

workflow PACBIO_BARCODE_CHECK {
    take:
    reference_tuple             // tuple    [[meta.id], reference ]
    pacbio_tuple                // tuple    [[meta.id], pacbio-files]
    barcodes                    // tuple    [[meta.id], barcode-file]
    barcode_multiplex           // val      (csv-list-string)

    main:
    ch_versions             = Channel.empty()

    //
    // MODULE: CHECK FOR KNOWN BARCODES IN SAMPLE DATA
    //
    CHECK_BARCODE (
        pacbio_tuple,
        barcodes.map{it[1]},
        barcode_multiplex
    )
    ch_versions     = ch_versions.mix(CHECK_BARCODE.out.versions)

    //
    // LOGIC: INCASE THE PIPELINE MANAGES TO CONTINUE AFTER FAILING CHECK_BARCODE
    //          HERE WE ENSURE THE REST OF THE SUBWORKFLOW DOES NOT RUN
    //
    CHECK_BARCODE.out.result
        .branch {
            valid   :   it.toString().contains('barcodes')
            invalid :   !it.toString().contains('barcodes')
        }
        .set { gatekeeping }

    //
    // LOGIC: ENSURE THE VALID CHANNEL IS MIXED WITH THE BARCODES CHANNEL
    //          ACTS AS A GATEKEEPER FOR THE FLOW
    //
    gatekeeping.valid
        .combine( barcodes )
        .map {str, meta, file ->
                file
        }
        .set {ch_new_barcodes}

    //
    // MODULE: GENERATE BLAST DB ON PACBIO BARCODES
    //
    BLAST_MAKEBLASTDB (
        ch_new_barcodes
    )
    ch_versions     = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

    //
    // MODULE: RUN BLAST WITH GENOME AGAINST BARCODE DB
    //
    BLAST_BLASTN (
        reference_tuple,
        BLAST_MAKEBLASTDB.out.db
    )
    ch_versions     = ch_versions.mix(BLAST_BLASTN.out.versions)

    //
    // LOGIC: FOR I (MAPPED TO OTHER CHANNELS) IN CSV LIST RUN FILTER BLAST
    //
    barcode_multiplex
        .map { it ->
            tuple( it.split(',') )
        }
        .flatten()
        .combine( reference_tuple )
        .combine( BLAST_BLASTN.out.txt )
        .multiMap { code, ref_meta, ref, blast_meta, blast ->
            barcodes: code
            reference: tuple( ref_meta, ref )
            blastdata: tuple( blast_meta, blast )
        }
        .set {testing}

    //
    // MODULE: CREATE A FILTERED BLAST OUTPUT PER BARCODE
    //
    FILTER_BARCODE (
        testing.reference,
        testing.blastdata,
        testing.barcodes
    )
    ch_versions     = ch_versions.mix(FILTER_BARCODE.out.versions)

    emit:
    filtered        = FILTER_BARCODE.out.debarcoded
    versions        = ch_versions.ifEmpty(null)
}
