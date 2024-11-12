include { TRAILINGNS      } from '../../modules/local/trailingns'

workflow TRAILINGNS_CHECK {

    take:
    reference_tuple             // tuple [ val(meta), path(fasta) ]

    main:

    ch_versions         = Channel.empty()

    //
    // MODULE: TRIM LENGTHS OF N'S FROM THE INPUT GENOME AND GENERATE A REPORT ON LENGTH AND LOCATION
    //
    TRAILINGNS( reference_tuple )
    ch_versions         = ch_versions.mix( TRAILINGNS.out.versions )

    emit:
    trailing_ns_report  = TRAILINGNS.out.trailing_ns_report
    versions            = ch_versions.ifEmpty( null )
}
