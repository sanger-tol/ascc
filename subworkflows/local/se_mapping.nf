include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_SE        } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_MERGE                             } from '../../modules/nf-core/samtools/merge/main'

workflow SE_MAPPING {

    take:
    reference_data_tuple     // Channel [ val(meta), path(file), path(file) ]

    main:
    ch_versions     = Channel.empty()
    ch_align_bams   = Channel.empty()

    //
    // LOGIC: MAKE MINIMAP INPUT CHANNEL AND MAKE BRANCHES BASED ON INPUT READ TYPE
    //
    reference_data_tuple
        .map { meta, ref, reads_path ->
            tuple(
                [   id          : meta.id,
                    single_end  : true,
                    readtype    : params.reads_type
                ],
                ref,
                reads_path,
                "bai",
                true,
                false,
                false,
                params.reads_type
            )
        }
        .multiMap{
            meta, reference, reads, index_format, bam_output, cigar_paf, cigar_bam, read_type ->
            reads_ch: tuple(meta, reads)
            refer_ch: tuple(meta, reference)
            inx_frmt: index_format
            bam_outp: bam_output
            cigar_pf: cigar_paf
            cigar_bm: cigar_bam
        }
        .set { minimap_se_input }


    //
    // MODULE: SINGLE END MAPPING THE READS WITH MINIMAP2
    //
    MINIMAP2_ALIGN_SE (
            minimap_se_input.reads_ch,
            minimap_se_input.refer_ch,
            minimap_se_input.bam_outp,
            minimap_se_input.inx_frmt,
            minimap_se_input.cigar_pf,
            minimap_se_input.cigar_bm
    )


    //
    // LOGIC: COLLECT THE OUTPUT FROM MINIMAP2
    //
    MINIMAP2_ALIGN_SE.out.bam
        .groupTuple(by: 0)
        .map{ meta, files ->
            tuple( meta, files.flatten())
        }
        .set { collected_files_for_merge }


    //
    // MODULE: MERGE ALL OUTPUT BAM
    //
    SAMTOOLS_MERGE(
        collected_files_for_merge,
        [[],[]],
        [[],[]]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)


    emit:
    versions       = ch_versions
    mapped_bam     = SAMTOOLS_MERGE.out.bam
}
