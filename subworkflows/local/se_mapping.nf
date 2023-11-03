include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_HIFI        } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_CLR         } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_ONT         } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_MERGE                                } from '../../modules/nf-core/samtools/merge/main'

workflow SE_MAPPING {

    take:
    reference_tuple          // Channel [ val(meta), path(file) ]
    assembly_path            // Channel path(file)
    pacbio_tuple             // Channel [ val(meta), val( str ) ]
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
                    split_prefix: read_path.toString().split('/')[-1].split('.fa.gz')[0]
                ],
                read_path,
                ref,
                true,
                false,
                false,
                platform
            )
        }
        .branch {
            hifi               : it[6] == "hifi"
            clr                : it[6] == "clr"
            ont                : it[6] == "ont"
            }
        .set { minimap_se_input }

    //
    // PROCESS: MULTIMAP TO MAKE BOOLEAN ARGUMENTS FOR MINIMAP HIFI MAPPING INPUT
    //
    minimap_se_input.hifi
        .multiMap { meta, read_path, ref, bam_output, cigar_paf, cigar_bam, platform ->
            read_tuple          : tuple( meta, read_path)
            ref                 : ref
            bool_bam_ouput      : bam_output
            bool_cigar_paf      : cigar_paf
            bool_cigar_bam      : cigar_bam
        }
        .set { hifi }

    //
    // PROCESS: MULTIMAP TO MAKE BOOLEAN ARGUMENTS FOR MINIMAP CLR MAPPING INPUT
    //
    minimap_se_input.clr
        .multiMap { meta, read_path, ref, bam_output, cigar_paf, cigar_bam, platform ->
            read_tuple          : tuple( meta, read_path)
            ref                 : ref
            bool_bam_ouput      : bam_output
            bool_cigar_paf      : cigar_paf
            bool_cigar_bam      : cigar_bam
        }
        .set { clr }

    //
    // PROCESS: MULTIMAP TO MAKE BOOLEAN ARGUMENTS FOR MINIMAP ONT MAPPING INPUT
    //
    minimap_se_input.ont
        .multiMap { meta, read_path, ref, bam_output, cigar_paf, cigar_bam, platform ->
            read_tuple          : tuple( meta, read_path)
            ref                 : ref
            bool_bam_ouput      : bam_output
            bool_cigar_paf      : cigar_paf
            bool_cigar_bam      : cigar_bam
        }
        .set { ont }

    //
    // MOUDLES: MAPPING DIFFERENT TYPE OF READ AGAINIST REFERENCE
    //
    if ( platform.filter { it == "hifi" } ){
        MINIMAP2_ALIGN_HIFI (
            hifi.read_tuple,
            hifi.ref,
            hifi.bool_bam_ouput,
            hifi.bool_cigar_paf,
            hifi.bool_cigar_bam
        )
        ch_bams = MINIMAP2_ALIGN_HIFI.out.bam
    } 
    else if ( platform.filter { it == "clr" } ){
        MINIMAP2_ALIGN_CLR (
            clr.read_tuple,
            clr.ref,
            clr.bool_bam_ouput,
            clr.bool_cigar_paf,
            clr.bool_cigar_bam
        )
        ch_bams = MINIMAP2_ALIGN_CLR.out.bam
    }
    else if ( platform.filter { it == "ont" } ){
        MINIMAP2_ALIGN_ONT (
            ont.read_tuple,
            ont.ref,
            ont.bool_bam_ouput,
            ont.bool_cigar_paf,
            ont.bool_cigar_bam
        )
        ch_bams = MINIMAP2_ALIGN_ONT.out.bam
    }

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
    tuple val(meta), path("in/*.{fa,fasta}.{gz}")

    "true"
}