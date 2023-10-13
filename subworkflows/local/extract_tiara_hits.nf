include { TIARA_TIARA } from '../../modules/nf-core/tiara/tiara/main'

workflow EXTRACT_TIARA_HITS {

    take:
    reference_tuple     // Channel [ val(meta), path(file) ]

    main:
    ch_versions     = Channel.empty()

    TIARA_TIARA (
        reference_tuple
    )
    ch_versions     = ch_versions.mix( TIARA_TIARA.out.versions )

    emit:
    ch_tiara            = TIARA_TIARA.out.classifications
    versions            = ch_versions.ifEmpty(null)

}
