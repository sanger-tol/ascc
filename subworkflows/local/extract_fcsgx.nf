include { FCS_FCSGX } from '../../modules/nf-core/fcs/fcsgx/main'



workflow EXTRACT_FCSGX {

    take:
    reference_tuple     // Channel [ val(meta), path(file) ]
    fcsgxpath

    main:
    ch_versions     = Channel.empty()

    Channel
    .of('all.gxi', 'all.gxs',  'all.taxa.tsv', 'all.meta.jsonl', 'all.blast_div.tsv.gz')
    .combine(fcsgxpath)
    .map {suxfix, dbpath -> [file(dbpath + '/' + suxfix)]}
    .collect()
    .set {fcsgxdb}

    fcsgxdb.view()

    FCS_FCSGX (
        reference_tuple,
        fcsgxdb
    )
    ch_versions     = ch_versions.mix( FCS_FCSGX.out.versions )

    emit:
    ch_fcs_gx_report    = FCS_FCSGX.out.fcs_gx_report
    ch_taxonomy_report  = FCS_FCSGX.out.taxonomy_report
    versions            = ch_versions.ifEmpty(null)

}