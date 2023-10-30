include { SE_MAPPING                                    } from './se_mapping'
include { PE_MAPPING as PE_MAPPING_ILLUMINA             } from './pe_mapping'
include { SAMTOOLS_MERGE                                } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX                                } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_DEPTH                                } from '../../modules/nf-core/samtools/depth/main'

workflow RUN_READ_COVERAGE {

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

    SE_MAPPING (
        reference_tuple,
        assembly_path,
        pacbio_tuple,
        platform
    )
    ch_versions = ch_versions.mix(SE_MAPPING.out.versions)
    ch_align_bams = SE_MAPPING.out.bam

    PE_MAPPING_ILLUMINA  (
        reference_tuple,
        assembly_path,
        pacbio_tuple,
        platform
    )
    ch_versions = ch_versions.mix(PE_MAPPING.out.versions)

    ch_align_bams
        .mix( PE_MAPPING.out.bam )
        .set { ch_bams }

    
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

    SAMTOOLS_MERGE(
        collected_files_for_merge,
        reference_tuple,
        [[],[]]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)


    SAMTOOLS_SORT (
        SAMTOOLS_MERGE.out.bam
    )
    ch_versions = ch_versions.mix( SAMTOOLS_MERGE.out.versions )

    SAMTOOLS_INDEX (
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix( SAMTOOLS_INDEX.out.versions )

    SAMTOOLS_DEPTH (
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix( SAMTOOLS_DEPTH.out.versions ) 


    emit:
    versions       = ch_versions.ifEmpty(null)
    //bam_ch         = MINIMAP2_ALIGN_HIFI.out.bam
}