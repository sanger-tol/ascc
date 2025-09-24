include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_ILLUMINA } from '../../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_MERGE                            } from '../../../modules/nf-core/samtools/merge/main'

workflow PE_MAPPING {

    take:
    reference_data_tuple     // Channel [ val(meta), path(file), path(file) ]

    main:
    ch_versions     = Channel.empty()

    //
    // LOGIC: MAKE MINIMAP INPUT CHANNEL
    //
    reference_data_tuple
        .map { meta, ref, reads_path, reads_type ->
            tuple(
                [   id          : meta.id,
                    single_end  : false,
                    readtype: reads_type.toString()
                ],
                reads_path,
                ref,
                true,
                false,
                false,
                reads_type
            )
        }
        .set { pe_input }


    //
    // LOGIC: MULTIMAP TO MAKE BOOLEAN ARGUMENTS
    //
    pe_input
        .multiMap { meta, reads_path, ref, bam_output, cigar_paf, cigar_bam, reads_type ->
            read_tuple          : tuple( meta, read_path)
            ref                 : tuple( meta, ref)
            bam_index_extension : "csi"
            bool_bam_ouput      : bam_output
            bool_cigar_paf      : cigar_paf
            bool_cigar_bam      : cigar_bam
        }
        .set { illumina_input }


    //
    // MODULE: PAIRED END READ MAPPING USING MINIMAP
    //
    MINIMAP2_ALIGN_ILLUMINA (
        illumina_input.read_tuple,
        illumina_input.ref,
        illumina_input.bool_bam_ouput,
        illumina_input.bam_index_extension,
        illumina_input.bool_cigar_paf,
        illumina_input.bool_cigar_bam
    )
    ch_versions     = ch_versions.mix(MINIMAP2_ALIGN_ILLUMINA.out.versions)


    MINIMAP2_ALIGN_ILLUMINA.out.bam
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
        reference_data_tuple,
        [[],[]],
        [[],[]]
    )
    ch_versions     = ch_versions.mix(SAMTOOLS_MERGE.out.versions)


    emit:
    versions        = ch_versions
    mapped_bam      = SAMTOOLS_MERGE.out.bam
}
