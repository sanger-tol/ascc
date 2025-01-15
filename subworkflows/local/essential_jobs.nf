
include { FILTER_FASTA                                  } from '../../modules/local/filter_fasta'
include { GC_CONTENT                                    } from '../../modules/local/gc_content'
include { GUNZIP                                        } from '../../modules/nf-core/gunzip/main'

include { GENERATE_GENOME                               } from '../../subworkflows/local/generate_genome'
include { TRAILINGNS_CHECK                              } from '../../subworkflows/local/trailingns_check'

workflow ESSENTIAL_JOBS {

    take:
    input_ref   // Channel [ val(meta), path(file) ]

    main:
    ch_versions             = Channel.empty()


    //
    // LOGIC: GUNZIP INPUT DATA IF GZIPPED, OTHERWISE PASS
    //
    input_ref
        .map { meta, ref ->
            def unzipped
            if (ref.name.endsWith('.gz')) {
                processed_file = ref.uncompress()
                println "Unzipped: ${ref.name} to ${processed_file.name}"
            } else {
                processed_file = ref
                println "No unzipping needed for: ${ref.name}"
            }
            [ meta, file(processed_file) ]
        }
        .set {ch_unzipped}

    println "Things ${ch_unzipped}"
    //
    // LOGIC: INJECT SLIDING WINDOW VALUES INTO REFERENCE
    //
    ch_unzipped
        .map { meta, ref ->
            tuple([ id      : meta.id,
                    sliding : params.seqkit_sliding,
                    window  : params.seqkit_window,
                    taxid   : params.taxid
                ],
                ref
            )
        }
        .set { new_input_fasta }

    println "Injected ${new_input_fasta}"


    //
    // MODULE: FILTER THE INPUT FASTA FOR LENGTHS OF SEQUENCE BELOW A CONFIG VARIABLE
    //
    FILTER_FASTA(
        new_input_fasta
    )
    ch_versions             = ch_versions.mix(FILTER_FASTA.out.versions)


    //
    // MODULE: CALCULATE GC CONTENT PER SCAFFOLD IN INPUT FASTA
    //
    GC_CONTENT (
        FILTER_FASTA.out.fasta
    )
    ch_versions             = ch_versions.mix(GC_CONTENT.out.versions)


    //
    // SUBWORKFLOW: GENERATE GENOME FILE - NA
    //
    GENERATE_GENOME (
        FILTER_FASTA.out.fasta,
        params.pacbio_barcode_names
    )
    ch_versions             = ch_versions.mix(GENERATE_GENOME.out.versions)


    //
    // SUBWORKFLOW: GENERATE A REPORT ON LENGTHS OF N's IN THE INPUT GENOME
    //
    TRAILINGNS_CHECK (
        FILTER_FASTA.out.fasta
    )
    ch_versions             = ch_versions.mix(TRAILINGNS_CHECK.out.versions)


    emit:
    reference_tuple_from_GG             = GENERATE_GENOME.out.reference_tuple
    reference_with_seqkit               = new_input_fasta
    dot_genome                          = GENERATE_GENOME.out.dot_genome
    gc_content_txt                      = GC_CONTENT.out.txt
    versions                            = ch_versions
}
