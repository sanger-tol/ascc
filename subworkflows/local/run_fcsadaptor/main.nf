//
// NF-CORE MODULE IMPORTS
//
include { FCS_FCSADAPTOR as FCS_FCSADAPTOR_PROK } from '../../../modules/nf-core/fcs/fcsadaptor/main'
include { FCS_FCSADAPTOR as FCS_FCSADAPTOR_EUK  } from '../../../modules/nf-core/fcs/fcsadaptor/main'

workflow RUN_FCSADAPTOR {
    take:
    reference_tuple

    main:
    ch_versions     = channel.empty()


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

    FCS_FCSADAPTOR_EUK.out.adaptor_report
        .map{ meta, file -> [meta.id, file] }
        .mix(
            FCS_FCSADAPTOR_PROK.out.adaptor_report
                .map{ meta, file -> [meta.id, file] }
        )
        .groupTuple()
        // NOTE: `groupTuple()` emits its internal ArrayBag, and that SAME instance is
        //       handed to every consumer of this channel. Any in-place operation
        //       (sort/unique/add) by one consumer corrupts it for the others, which
        //       surfaces as `Unexpected error [ConcurrentModificationException]`.
        //       Emit a defensive copy so each downstream branch owns its own list.
        .map { id, files -> [[id: id], new ArrayList(files)] }
        .set { ch_fcsadapt }

    emit:
    ch_joint_report = ch_fcsadapt
    ch_euk          = FCS_FCSADAPTOR_EUK.out.adaptor_report
    ch_prok         = FCS_FCSADAPTOR_PROK.out.adaptor_report

    versions        = ch_versions
}
