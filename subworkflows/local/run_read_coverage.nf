include { SE_MAPPING                                    } from './se_mapping'
include { PE_MAPPING as PE_MAPPING_ILLUMINA             } from './pe_mapping'
include { SAMTOOLS_MERGE                                } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX                                } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT                                 } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_DEPTH                                } from '../../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_DEPTH_AVERAGE_COVERAG                } from '../../modules/local/samtools_depth_average_coverage'

workflow RUN_READ_COVERAGE {

    take:
    reference_tuple          // Channel [ val(meta), path(file) ]
    assembly_path            // Channel path(file)
    pacbio_tuple             // Channel [ val(meta), val( str ) ]
    platform                // Channel val( str )

    main:
    ch_versions     = Channel.empty()
    ch_align_bam    = Channel.empty()

    
    //
    // MODULE: GETS PACBIO READ PATHS FROM READS_PATH
    //
    if ( platform.filter { it == "hifi" } || platform.filter { it == "clr" } || platform.filter { it == "ont" } ) { 
        SE_MAPPING (
            reference_tuple,
            assembly_path,
            pacbio_tuple,
            platform
        )
        ch_versions = ch_versions.mix(SE_MAPPING.out.versions)
        ch_align_bam
            .mix( SE_MAPPING.out.mapped_bam )
            .set { merged_bam }
    }
    else if ( platform.filter { it == "illumina" } ) { 

        PE_MAPPING_ILLUMINA  (
            reference_tuple,
            assembly_path,
            pacbio_tuple,
            platform
        )
        ch_versions = ch_versions.mix(PE_MAPPING_ILLUMINA.out.versions)
        ch_align_bam
            .mix( PE_MAPPING_ILLUMINA.out.mapped_bam )
            .set { merged_bam }
    }

    SAMTOOLS_SORT (
        merged_bam
    )
    ch_versions = ch_versions.mix( SAMTOOLS_SORT.out.versions )

    SAMTOOLS_INDEX (
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix( SAMTOOLS_INDEX.out.versions )

    SAMTOOLS_DEPTH (
        SAMTOOLS_SORT.out.bam,
        [[],[]]
    )
    ch_versions = ch_versions.mix( SAMTOOLS_DEPTH.out.versions )

    SAMTOOLS_DEPTH_AVERAGE_COVERAG (
        SAMTOOLS_DEPTH.out.tsv
    )
    ch_versions = ch_versions.mix( SAMTOOLS_DEPTH_AVERAGE_COVERAG.out.versions )

    emit:
    versions       = ch_versions.ifEmpty(null)
    tsv_ch         = SAMTOOLS_DEPTH_AVERAGE_COVERAG.out.average_coverage
}