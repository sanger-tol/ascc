include { SE_MAPPING                                    } from './se_mapping'
include { PE_MAPPING                                    } from './pe_mapping'
include { SAMTOOLS_MERGE                                } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX                                } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT                                 } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_DEPTH                                } from '../../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_DEPTH_AVERAGE_COVERAGE               } from '../../modules/local/samtools_depth_average_coverage'

workflow RUN_READ_COVERAGE {

    take:
    reference_tuple          // Channel [ val(meta), path(file) ]
    reads
    pacbio_data              // [Path(1), Path(2)...]
    platform                 // Channel val( str )

    main:
    ch_versions     = Channel.empty()
    ch_align_bam    = Channel.empty()
    ch_refer_bam    = Channel.empty()


    //
    // LOGIC: GETS PACBIO READ PATHS FROM READS_PATH
    //
    reads
        .flatten()
        .set{ collection_of_reads }

    reference_tuple
        .combine(collection_of_reads)
        .set { ref_and_data }

    //
    // LOGIC: CHECK IF THE INPUT READ FILE IS PAIRED END OR SINGLE END BASED ON THE READ PLATFORM
    // THEN RUN MINIMAP
    // - Removed the mix function from this as it is not needed, there shouldn't be multiple read
    // types
    //
    reference_tuple
        .map{meta, file ->
            tuple(meta, pacbio_data)
        }
        .set { pacbio_tuple }

    Channel
        .of(platform)
        .set {platform_type}

    if ( platform == "hifi" || platform == "clr" || platform == "ont" ) {

        //
        // MODULE: RUN SINGLE END MAPPING ON THE REFERENCE AND LONGREAD DATA
        //
        SE_MAPPING (
            ref_and_data
        )
        ch_versions = ch_versions.mix(SE_MAPPING.out.versions)

        SE_MAPPING.out.mapped_bam
            .set { ch_align_bam }

    }
    else if ( platform == "illumina" ) {

        //
        // MODULE: RUN PAIRED END MAPPING ON THE REFERENCE AND LONGREAD DATA
        //
        PE_MAPPING  (
            reference_tuple,
            pacbio_tuple,
            platform_type
        )
        ch_versions = ch_versions.mix(PE_MAPPING.out.versions)

        PE_MAPPING.out.mapped_bam
            .set { ch_align_bam }

    }


    //
    // MODULE: SORT THE MAPPED BAM
    //
    SAMTOOLS_SORT (
        ch_align_bam,
        [[],[]]
    )
    ch_versions = ch_versions.mix( SAMTOOLS_SORT.out.versions )


    //
    // MODULE: INDEX THE SORTED MAPPED BAM
    //
    SAMTOOLS_INDEX (
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix( SAMTOOLS_INDEX.out.versions )


    //
    // MODULE: GET READ DEPTH ACROSS THE GENOME
    //
    SAMTOOLS_DEPTH (
        SAMTOOLS_SORT.out.bam,
        [[],[]]
    )
    ch_versions = ch_versions.mix( SAMTOOLS_DEPTH.out.versions )


    //
    // MODULE: COMPUTE THE AVERAGE COVERAGE ACROSS EACH SCAFFOLD
    //
    SAMTOOLS_DEPTH_AVERAGE_COVERAGE (
        SAMTOOLS_DEPTH.out.tsv
    )
    ch_versions = ch_versions.mix( SAMTOOLS_DEPTH_AVERAGE_COVERAGE.out.versions )


    emit:
    tsv_ch          = SAMTOOLS_DEPTH_AVERAGE_COVERAGE.out.average_coverage
    bam_ch          = SAMTOOLS_SORT.out.bam
    versions        = ch_versions.ifEmpty(null)
}
