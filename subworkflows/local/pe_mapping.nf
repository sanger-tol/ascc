include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_ILLUMINA         } from '../../modules/nf-core/minimap2/align/main'

workflow PE_MAPPING {

    take:
    reference_tuple          // Channel [ val(meta), path(file) ]
    assembly_path            // Channel path(file)
    pacbio_tuple             // Channel [ val(meta), val( str ) ]
    platform                // Channel val( str )

    main:
    ch_versions     = Channel.empty()


    //
    // MODULE: GETS PACBIO READ PATHS FROM READS_PATH
    //
    ch_grabbed_read_paths       = GrabFiles( pacbio_tuple )

    ch_grabbed_read_paths
        .map { meta, files ->
            tuple( files )
        }
        .flatten()
        .set { ch_read_paths }

    reference_tuple
        .combine( ch_read_paths )
        .combine( read_type )
        .combine( se)
        .map { meta, ref, read_path, platform ->
            tuple(
                [   id          : meta.id,
                    single_end  : false,
                    split_prefix: read_path.toString().split('/')[-1].split('.fasta.gz')[0]
                ],
                read_path,
                ref,
                true,
                false,
                false,
                platform
            )
        }
        .set { pe_input }

    pe_input
        .multiMap { meta, read_path, ref, bam_output, cigar_paf, cigar_bam, platform ->
            read_tuple          : tuple( meta, read_path)
            ref                 : ref
            bool_bam_ouput      : bam_output
            bool_cigar_paf      : cigar_paf
            bool_cigar_bam      : cigar_bam
        }
        .set { illumina_input }

    MINIMAP2_ALIGN_ILLUMINA (
        illumina_input.read_tuple,
        illumina_input.ref,
        illumina_input.bool_bam_ouput,
        illumina_input.bool_cigar_paf,
        illumina_input.bool_cigar_bam
    )

    emit:
    versions           = ch_versions.ifEmpty(null)
    illumina_bam_ch    = MINIMAP2_ALIGN_ILLUMINA.out.bam
}

process GrabFiles {
    tag "${meta.id}"
    executor 'local'

    input:
    tuple val(meta), path("in")

    output:
    tuple val(meta), path("in/*fasta.gz")

    "true"
}