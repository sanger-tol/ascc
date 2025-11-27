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


    //
    // LOGIC: AT THIS POINT THE META CONTAINS JUNK THAT CAN 'CONTAMINATE' MATCHES,
    //          SO STRIP IT DOWN BEFORE USE, WE ALSO MERGE THE OUTPUT TOGETHER FOR SIMPLICITY
    //
    FCS_FCSADAPTOR_EUK.out.adaptor_report
        .combine(
            FCS_FCSADAPTOR_PROK.out.adaptor_report.map{meta, file -> file}
        )
        .map { meta, file1, file2 ->
            [[id: meta.id], file1, file2 ]
        }
        .set { ch_fcsadapt }

    emit:
    ch_joint_report = ch_fcsadapt
    ch_euk          = FCS_FCSADAPTOR_EUK.out.adaptor_report
    ch_prok         = FCS_FCSADAPTOR_PROK.out.adaptor_report

    versions        = ch_versions
}
