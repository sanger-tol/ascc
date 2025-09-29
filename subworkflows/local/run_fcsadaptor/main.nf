//
// NF-CORE MODULE IMPORTS
//
include { FCS_FCSADAPTOR as FCS_FCSADAPTOR_PROK } from '../../../modules/nf-core/fcs/fcsadaptor/main'
include { FCS_FCSADAPTOR as FCS_FCSADAPTOR_EUK  } from '../../../modules/nf-core/fcs/fcsadaptor/main'

workflow RUN_FCSADAPTOR {
    take:
    reference_tuple

    main:
    ch_versions     = Channel.empty()


    //
    // MODULE: FCS_FCSADAPTOR run on assembly fasta for prokaryote.
    //
    FCS_FCSADAPTOR_PROK (
        reference_tuple      // val(meta), path(fasta)
    )
    ch_versions     = ch_versions.mix(FCS_FCSADAPTOR_PROK.out.versions)


    //
    // MODULE: FCS_FCSADAPTOR run with eukaryote settings
    //
    FCS_FCSADAPTOR_EUK (
        reference_tuple      // val(meta), path(fasta)
    )
    ch_versions     = ch_versions.mix(FCS_FCSADAPTOR_EUK.out.versions)

    emit:
    ch_euk          = FCS_FCSADAPTOR_EUK.out.adaptor_report
    ch_prok         = FCS_FCSADAPTOR_PROK.out.adaptor_report

    versions        = ch_versions
}
