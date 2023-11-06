include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_SE        } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_MERGE                             } from '../../modules/nf-core/samtools/merge/main'

workflow SE_MAPPING {

    take:
    reference_tuple          // Channel [ val(meta), path(file) ]
    assembly_path            // Channel path(file)
    pacbio_tuple             // Channel [ val(meta), path(file) ]
    platform                // Channel val( str )

    main:
    ch_versions     = Channel.empty()
    ch_align_bams   = Channel.empty()

    //
    // PROCESS: GETS PACBIO READ PATHS FROM READS_PATH
    //
    ch_grabbed_read_paths       = GrabFiles( pacbio_tuple )

    ch_grabbed_read_paths
        .map { meta, files ->
            tuple( files )
        }
        .flatten()
        .set { ch_read_paths }

    //
    // PROCESS: MAKE MINIMAP INPUT CHANNEL AND MAKE BRANCHES BASED ON INPUT READ TYPE
    //
    reference_tuple
        .combine( ch_read_paths )
        .combine( platform )
        .map { meta, ref, read_path, platform ->
            tuple(
                [   id          : meta.id,
                    single_end  : true,
                    readtype    : platform.toString()
                ],
                read_path,
                ref,
                true,
                false,
                false,
                platform
            )
        }
        .set { minimap_se_input }

    //
    // PROCESS: MULTIMAP TO MAKE BOOLEAN ARGUMENTS FOR MINIMAP HIFI MAPPING INPUT
    //
    minimap_se_input
        .multiMap { meta, read_path, ref, bam_output, cigar_paf, cigar_bam, platform ->
            read_tuple          : tuple( meta, read_path)
            ref                 : ref
            bool_bam_ouput      : bam_output
            bool_cigar_paf      : cigar_paf
            bool_cigar_bam      : cigar_bam
        }
        .set { se_input }

    //
    // MOUDLES: MAPPING DIFFERENT TYPE OF READ AGAINIST REFERENCE
    //

    MINIMAP2_ALIGN_SE (
            se_input.read_tuple,
            se_input.ref,
            se_input.bool_bam_ouput,
            se_input.bool_cigar_paf,
            se_input.bool_cigar_bam
    )
    ch_bams = MINIMAP2_ALIGN_SE.out.bam
   

    ch_bams
        .map { meta, file ->
            tuple( file )
        }
        .collect()
        .map { file ->
            tuple (
                [ id    : file[0].toString().split('/')[-1].split('_')[0] ], // Change sample ID
                file
            )
        }
        .set { collected_files_for_merge }

    //
    // MODULE: MERGE ALL OUTPUT BAM
    //
    SAMTOOLS_MERGE(
        collected_files_for_merge,
        reference_tuple,
        [[],[]]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    emit:
    versions       = ch_versions.ifEmpty(null)
    mapped_bam     = SAMTOOLS_MERGE.out.bam
}

process GrabFiles {
    tag "${meta.id}"
    executor 'local'

    input:
    tuple val(meta), path("in")

    output:
    tuple val(meta), path("in/*.fa.gz")

    "true"
}