// LOCAL IMPORTS
include { FILTER_FASTA   } from '../../../modules/local/filter/fasta/main'
include { GC_CONTENT     } from '../../../modules/local/gc/content/main'
include { TRAILINGNS     } from '../../../modules/local/trailingns/trailingns/main'

// NF-CORE IMPORTS
include { GNU_SORT       } from '../../../modules/nf-core/gnu/sort/main'
include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx/main'


workflow ESSENTIAL_JOBS {
    take:
    input_ref // channel [ val(meta), path(file) ]

    main:
    ch_versions = channel.empty()
    //
    // MODULE: FILTER/BREAK THE INPUT FASTA FOR LENGTHS OF SEQUENCE BELOW A 1.9Gb THRESHOLD, MORE THAN THIS WILL BREAK SOME TOOLS
    //
    // Determine if FCS-adaptor will be run based on run_fcs_adaptor parameter
    def run_fcs_adaptor = (params.run_fcs_adaptor == "both" || (params.run_fcs_adaptor == "genomic" && params.genomic_only) || (params.run_fcs_adaptor == "organellar" && !params.genomic_only))

    FILTER_FASTA(
        input_ref,
        run_fcs_adaptor,
    )
    filter_fasta_sanitation_log         = FILTER_FASTA.out.sanitation_log
                                             .map{ meta, _file -> tuple([id: meta.id], _file) }
    filter_fasta_length_filtering_log   = FILTER_FASTA.out.length_filtering_log
                                             .map{ meta, _file -> tuple([id: meta.id], _file) }

    //
    // MODULE: CALCULATE GC CONTENT PER SCAFFOLD IN INPUT FASTA
    //
    GC_CONTENT(
        FILTER_FASTA.out.fasta
    )


    //
    // MODULE: GENERATE INDEX OF REFERENCE
    //          EMITS REFERENCE INDEX FILE MODIFIED FOR SCAFF SIZES
    //
    SAMTOOLS_FAIDX(
        FILTER_FASTA.out.fasta.map { meta, fasta -> tuple(meta, fasta, []) },
        true,
    )


    //
    // MODULE: SORT CHROM SIZES BY CHOM SIZE NOT NAME
    //
    GNU_SORT(
        SAMTOOLS_FAIDX.out.sizes.map { meta, _file -> tuple(meta, _file, "sizes") }
    )


    //
    // MODULE: TRIM LENGTHS OF N'S FROM THE INPUT GENOME AND GENERATE A REPORT ON LENGTH
    //          AND LOCATION
    //
    TRAILINGNS(
        FILTER_FASTA.out.fasta
    )
    trailing_ns_report  = TRAILINGNS.out.trailing_ns_report
                                .map { meta, _file -> tuple([id: meta.id], _file) }

    emit:
    reference_tuple                   = FILTER_FASTA.out.fasta
    dot_genome                        = GNU_SORT.out.sorted
    gc_content_txt                    = GC_CONTENT.out.txt
    trailing_ns_report
    filter_fasta_sanitation_log
    filter_fasta_length_filtering_log
}
