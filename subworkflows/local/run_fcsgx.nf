include { FCS_FCSGX             } from '../../modules/nf-core/fcs/fcsgx/main'
include { PARSE_FCSGX_RESULT    } from '../../modules/local/parse_fcsgx_result'

workflow RUN_FCSGX {

    take:
    reference               // Channel [ val(meta), path(file) ]
    fcsgxpath               // Channel path(file)
    taxid                   // Channel val(taxid)
    ncbi_rankedlineage_path // Channel path(file)

    main:
    ch_versions     = Channel.empty()


    //
    // LOGIC: Create input channel for FCS_FCSGX, taxid is required to be the meta id.
    //
    reference
        .combine(
            Channel.of(taxid)
        )
        .map { it ->
                tuple ([    id:     it[0].id,
                            taxid:  it[2]       ],
                        it[1])
            }
        .set { reference_with_taxid }


    //
    // MODULE: FCS_FCSGX RUN ON ASSEMBLY FASTA TUPLE WITH THE TAXID AGAINST THE FCSGXDB
    //
    FCS_FCSGX (
        reference_with_taxid,
        fcsgxpath
    )
    ch_versions     = ch_versions.mix( FCS_FCSGX.out.versions )


    //
    // MODULE: CREATE INPUT CHANNEL FOR PARSING RESULT MODULE
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
