include { SE_MAPPING                } from './se_mapping'
//include { SE_MAPPING as SE_MAPPING_CLR                  } from './se_mapping'
//include { SE_MAPPING as SE_MAPPING_ONT                  } from './se_mapping'
//include { PE_MAPPING as PE_MAPPING_ILLUMINA             } from './pe_mapping'

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

    PE_MAPPING (
        reference_tuple,
        assembly_path,
        pacbio_tuple,
        platform
    )


    emit:
    versions       = ch_versions.ifEmpty(null)
    //bam_ch         = MINIMAP2_ALIGN_HIFI.out.bam
}