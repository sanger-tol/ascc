//
// PACBIO_BARCODE_CHECK IDENTIFIED LOCATIONS OF BARCODE SEQUENCES IN THE INPUT ASSEMBLY
//

//
// MODULE IMPORT BLOCK
//
include { CHECK_BARCODE                                         } from '../../../modules/local/check/barcode/main'
//include { BLAST_MAKEBLASTDB_BARCODES as BLAST_MAKEBLASTDB       } from '../../../modules/local/makeblastdb_barcodes'
include { BLAST_MAKEBLASTDB                                     } from '../../../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN                                          } from '../../../modules/nf-core/blast/blastn/main'
include { FILTER_BARCODE                                        } from '../../../modules/local/filter/barcode/main'

workflow PACBIO_BARCODE_CHECK {
    take:
    reference_tuple         // tuple    [[meta.id], reference ]
    pacbio_data             // tuple    pacbio-files-folder
    pacbio_type             // val      (params.pacbio_type)
    barcodes_file           // tuple    [[meta.id], barcode-file]
    barcode_names           // val      (csv-list-string)

    main:
    ch_versions             = Channel.empty()

    //
    // MODULE: CHECK FOR KNOWN BARCODES IN SAMPLE DATA
    //
    reference_tuple
        .map{
            tuple([id: it[0].id], pacbio_data)
        }
        .set {pacbio_data_tuple}

    CHECK_BARCODE (
        pacbio_data_tuple,
        barcodes_file,
        barcode_names
    )
    ch_versions     = ch_versions.mix(CHECK_BARCODE.out.versions)


    //
    // LOGIC: ENSURE THE VALID CHANNEL IS MIXED WITH THE BARCODES CHANNEL
    //          ACTS AS A GATEKEEPER FOR THE FLOW
    //
    Channel.fromPath(barcodes_file)
        .set {barcode_data_file}

    CHECK_BARCODE.out.result
        .filter{it.contains("barcodes")} // Indicates it is a valid barcode
        .combine( barcode_data_file )
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
