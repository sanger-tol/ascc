//
// PACBIO_BARCODE_CHECK IDENTIFIED LOCATIONS OF BARCODE SEQUENCES IN THE INPUT ASSEMBLY
//

//
// MODULE IMPORT BLOCK
//
include { CHECK_BARCODE          } from '../../modules/local/check_barcode'
include { BLAST_BLASTN           } from '../../modules/nf-core/blast/blastn'
include { FILTER_BARCODE         } from '../../modules/local/filter_barcode'

workflow PACBIO_BARCODE_CHECK {
    take:
    reference_tuple             // tuple    [[meta.id], reference ]
    pacbio_tuple                // tuple    [[meta.id], pacbio-files]
    barcode_multiplex           // val      (csv-list-string)
    blast_db                    // tuple    [[meta.id], blast_db]
    barcode_file

    main:
    ch_versions             = Channel.empty()

    //
    // MODULE: CHECK FOR KNOWN BARCODES IN SAMPLE DATA
    //

    CHECK_BARCODE (
        pacbio_tuple,
        barcode_file,
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
        .combine( barcode_file )
        .map {str, meta, file ->
                file
        }
        .set {ch_new_barcodes}

    //
    // MODULE: RUN BLAST WITH GENOME AGAINST BARCODE DB
    //

    reference_tuple
        .groupTuple(by: [1])
        .map{ it -> tuple( it[0][0], it[-1])} // get meta_1 and last item which will be the fasta
        .combine(blast_db)
        .multiMap{ meta, ref, meta_2, blast_db ->
            ref: tuple(meta, ref)
            blast: tuple(meta_2, blast_db)
        }
        .set{ done }
    
    // NOTE: WE ONLY NEED THIS TO RUN 1 PER INPUT GENOME
    BLAST_BLASTN (
        done.ref,
        done.blast
    )
    ch_versions     = ch_versions.mix(BLAST_BLASTN.out.versions)

    done.ref.view{"REFIES: $it"}

    //
    // MODULE: CREATE A FILTERED BLAST OUTPUT PER BARCODE
    //          MODULE RUNS ONCE PER REFERENCE AND INTERATES OVER BARCODE INTERNALLY
    //

    //NOTE: THE BELOW MULTIPLIED THE BARCODES BY BLAST_OUTPUT
    BLAST_BLASTN.out.txt
        .combine(barcode_multiplex)
        .multiMap{ meta, blast, codes ->
            blast: tuple(meta, blast)
            barcods: codes
        }
        .set {things}

    FILTER_BARCODE (
        done.ref,
        things.blast,
        things.barcods
    )
    ch_versions     = ch_versions.mix(FILTER_BARCODE.out.versions)

    FILTER_BARCODE.out.debarcoded.view{"OUTPUT: $it"}

    emit:
    // filtered        = FILTER_BARCODE.out.debarcoded
    versions        = ch_versions.ifEmpty(null)
}
