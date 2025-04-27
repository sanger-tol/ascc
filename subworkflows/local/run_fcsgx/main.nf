include { FCSGX_RUNGX           } from '../../../modules/nf-core/fcsgx/rungx/main'
include { PARSE_FCSGX_RESULT    } from '../../../modules/local/fcsgx/parse_results/main'

workflow RUN_FCSGX {

    take:
    reference               // Channel [ val(meta), path(file) ]
    fcsgxpath               // Channel path(file)
    ncbi_rankedlineage_path // Channel path(file)

    main:
    ch_versions     = Channel.empty()


    //
    // MODULE: FCSGX_RUNGX RUN ON ASSEMBLY FASTA TUPLE WITH THE TAXID AGAINST THE FCSGXDB
    //
    FCSGX_RUNGX (
        reference,
        fcsgxpath,
        []
    )
    ch_versions     = ch_versions.mix( FCSGX_RUNGX.out.versions )


    //
    // MODULE: CREATE INPUT CHANNEL FOR PARSING RESULT MODULE
    //
    FCSGX_RUNGX.out.fcsgx_report
        .map{ it ->
                tuple(  it[0],
                        it[1].getParent()
                )
        }
        .set { report_path }


    //
    // MODULE: PARSE_FCSGX_RESULT to parse the FCSGX_RUNGX result output in csv format.
    //
    PARSE_FCSGX_RESULT (
        report_path,
        ncbi_rankedlineage_path
    )
    ch_versions     = ch_versions.mix( PARSE_FCSGX_RESULT.out.versions )

    emit:
    fcsgxresult    = PARSE_FCSGX_RESULT.out.fcsgxresult
    fcsgx_report   = FCSGX_RUNGX.out.fcsgx_report
    taxonomy_report = FCSGX_RUNGX.out.taxonomy_report
    versions       = ch_versions.ifEmpty(null)

}
