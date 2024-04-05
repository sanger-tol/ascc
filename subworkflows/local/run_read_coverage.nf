include { SE_MAPPING                                    } from './se_mapping'
include { PE_MAPPING                                    } from './pe_mapping'
include { SAMTOOLS_MERGE                                } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX                                } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT                                 } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_DEPTH                                } from '../../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_DEPTH_AVERAGE_COVERAGE                } from '../../modules/local/samtools_depth_average_coverage'

workflow RUN_READ_COVERAGE {

    take:
    reference_tuple          // Channel [ val(meta), path(file) ]
    assembly_path            // Channel path(file)
    pacbio_tuple             // Channel [ val(meta), val( str ) ]
    platform                 // Channel val( str )

    main:
    ch_versions     = Channel.empty()
    ch_align_bam    = Channel.empty()


    //
    // LOGIC: CHECK IF THE INPUT READ FILE IS PAIRED END OR SINGLE END BASED ON THE READ PLATFORM, THEN RUN MINIMAP
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

        PE_MAPPING  (
            reference_tuple,
            assembly_path,
            pacbio_tuple,
            platform
        )
        ch_versions = ch_versions.mix(PE_MAPPING.out.versions)
        ch_align_bam
            .mix( PE_MAPPING.out.mapped_bam )
            .set { merged_bam }
    }

    //
    // MODULE: SORT MAPPED BAM
    //
    SAMTOOLS_SORT (
        merged_bam
    )
    ch_versions = ch_versions.mix( SAMTOOLS_SORT.out.versions )

    //
    // MODULE: INDEXING SORTED MAPPED BAM
    //
    SAMTOOLS_INDEX (
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix( SAMTOOLS_INDEX.out.versions )

    //
    // MODULE: GET READ DEPTH
    //
    SAMTOOLS_DEPTH (
        SAMTOOLS_SORT.out.bam,
        [[],[]]
    )
    ch_versions = ch_versions.mix( SAMTOOLS_DEPTH.out.versions )


    //
    // MODULE: COMPUTE THE AVERAGE COVERAGE
    //
    SAMTOOLS_DEPTH_AVERAGE_COVERAGE (
        SAMTOOLS_DEPTH.out.tsv
    )
    ch_versions = ch_versions.mix( SAMTOOLS_DEPTH_AVERAGE_COVERAGE.out.versions )

    emit:
    tsv_ch         = SAMTOOLS_DEPTH_AVERAGE_COVERAGE.out.average_coverage
    versions       = ch_versions.ifEmpty(null)
}