include { SAMTOOLS_DICT         } from '../../../modules/nf-core/samtools/dict/main'
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
    // MODULE: Use SAMTOOLS_DICT to get origin file and md5sum of each sequence
    //
    SAMTOOLS_DICT(
        reference
    )
    ch_versions     = ch_versions.mix( SAMTOOLS_DICT.out.versions )


    //
    // MODULE: FCSGX_RUNGX RUN ON ASSEMBLY FASTA TUPLE WITH THE TAXID AGAINST THE FCSGXDB
    //         PRODUCTION_MODE WILL CHANGE HOW THE FCSGX MODULE IS RUN E.G IT IS SPECIFIC FOR `module`
    //
    SAMTOOLS_DICT.out.dict
        .map { meta, ref, dict ->
            tuple(meta, ref)
        }
        .set { samtools_reference }

    FCSGX_RUNGX (
        samtools_reference,
        fcsgxpath,
        [],
        params.production_mode
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
    genomedict     = samtools_reference
    versions       = ch_versions

}
