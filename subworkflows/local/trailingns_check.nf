include { TRAILINGNS      } from '../../modules/local/trailingns'
//TRAILINGNS

workflow TRAILINGNS_CHECK {

    take:
    reference_tuple             // val(meta), path(fasta)

    main:

    ch_versions = Channel.empty()

    TRAILINGNS( reference_tuple )
    ch_versions = ch_versions.mix( TRAILINGNS.out.versions )

    //versions = ch_versions                     // channel: [ versions.yml ]
    emit:
    trailing_ns_report     = TRAILINGNS.out.trailing_ns_report
    versions                    = ch_versions.ifEmpty( null ) // channel: [ versions.yml ]
}

