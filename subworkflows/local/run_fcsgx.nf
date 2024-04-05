include { FCS_FCSGX }          from '../../modules/nf-core/fcs/fcsgx/main'
include { PARSE_FCSGX_RESULT } from '../../modules/local/parse_fcsgx_result'

workflow RUN_FCSGX {

    take:
    reference               // Channel [ val(meta), path(file) ]
    fcsgxpath               // Channel path(file)
    taxid                   // Channel val(taxid)
    ncbi_rankedlineage_path // Channel path(file)

    main:
    ch_versions     = Channel.empty()

    Channel
        .of(
            'all.gxi', 'all.gxs',  'all.taxa.tsv', 'all.meta.jsonl', 'all.blast_div.tsv.gz'
        )
        .combine(
            fcsgxpath
        )
        .map {suxfix, dbpath ->
            [file(dbpath + '/' + suxfix)]
        }
        .collect()
        .set { fcsgxdb }

    //
    // Create input channel for FCS_FCSGX, taxid is required to be the meta id.
    //
    reference
        .combine( taxid )
        .map { it ->
                tuple ([    id:     it[0].id,
                            taxid:  it[2]       ],
                        it[1])
            }
        .set { reference_with_taxid }

    //
    // MODULE: FCS_FCSGX run on assembly fasta tuple with taxid againist fcsgxdb.
    //
    FCS_FCSGX (
        reference_with_taxid,
        fcsgxdb
    )
    ch_versions     = ch_versions.mix( FCS_FCSGX.out.versions )

    //
    // Create input channel for parsing result module.
    //
    FCS_FCSGX.out.fcs_gx_report
        .map{ it ->
                tuple(  it[0],
                        it[1].getParent()
                )
        }
        .set { report_path }

    //
    // MODULE: PARSE_FCSGX_RESULT to parse the FCS_FCSGX result output in csv format.
    //
    PARSE_FCSGX_RESULT (
        report_path,
        ncbi_rankedlineage_path
    )
    ch_versions     = ch_versions.mix( PARSE_FCSGX_RESULT.out.versions )

    emit:
    fcsgxresult    = PARSE_FCSGX_RESULT.out.fcsgxresult
    versions       = ch_versions.ifEmpty(null)

}
