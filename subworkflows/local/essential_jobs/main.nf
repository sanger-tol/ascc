include { FILTER_FASTA                                  } from '../../../modules/local/filter/fasta/main'
include { GC_CONTENT                                    } from '../../../modules/local/gc/content/main'

include { GENERATE_GENOME                               } from '../../../subworkflows/local/generate_genomes/main'
include { TRAILINGNS_CHECK                              } from '../../../subworkflows/local/trailingns_check/main'

workflow ESSENTIAL_JOBS {

    take:
    input_ref   // Channel [ val(meta), path(file) ]

    main:
    ch_versions             = Channel.empty()


    //
    // LOGIC: INJECT SLIDING WINDOW VALUES INTO REFERENCE
    //
    input_ref
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


    //
    // MODULE: FILTER/BREAK THE INPUT FASTA FOR LENGTHS OF SEQUENCE BELOW A 1.9Gb THRESHOLD, MORE THAN THIS WILL BREAK SOME TOOLS
    //
    // Determine if FCS-adaptor will be run based on run_fcs_adaptor parameter
    def run_fcs_adaptor = (params.run_fcs_adaptor == "both" ||
                        (params.run_fcs_adaptor == "genomic" && params.genomic_only) ||
                        (params.run_fcs_adaptor == "organellar" && !params.genomic_only))

    FILTER_FASTA(
        new_input_fasta,
        run_fcs_adaptor
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
    trailing_ns_report                  = TRAILINGNS_CHECK.out.trailing_ns_report
    filter_fasta_sanitation_log         = FILTER_FASTA.out.sanitation_log.map { meta, file ->
        // Ensure we have a consistent structure
        [meta, file]
    }
    filter_fasta_length_filtering_log   = FILTER_FASTA.out.length_filtering_log.map { meta, file ->
        // Ensure we have a consistent structure
        [meta, file]
    }
    versions                            = ch_versions
}
