//
// MODULE IMPORT BLOCK
//
include { SAMTOOLS_DICT         } from '../../../modules/nf-core/samtools/dict/main'
include { FCSGX_RUNGX           } from '../../../modules/nf-core/fcsgx/rungx/main'
include { PARSE_FCSGX_RESULT    } from '../../../modules/local/fcsgx/parse_results/main'

workflow RUN_FCSGX {

    take:
    reference               // channel [ val(meta), path(file) ]
    fcsgxpath               // channel path(file)
    ncbi_rankedlineage_path // channel path(file)

    main:
    ch_versions     = channel.empty()

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
        .map { meta, ref, _dict ->
            [meta, ref]
        }
        .set { samtools_reference }

    FCSGX_RUNGX (
        samtools_reference,
        fcsgxpath,
        [],
        "production" in workflow.profile.tokenize(',')
    )
    ch_versions         = ch_versions.mix( FCSGX_RUNGX.out.versions )

    fcsgx_report_txt    = FCSGX_RUNGX.out.fcsgx_report
                            .map { meta, file ->
                                file ? [[ id: meta.id ], file] : [[:], []]
                            }

    fcsgx_taxonomy_rpt  = FCSGX_RUNGX.out.taxonomy_report
                            .map { meta, file ->
                                file ? [[ id: meta.id ], file] : [[:], []]
                            }

    //
    // MODULE: CREATE INPUT CHANNEL FOR PARSING RESULT MODULE
    //
    fcsgx_report_txt
        .map{ meta, file ->
            [meta, file.getParent()]
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

    fcsgxresult     = PARSE_FCSGX_RESULT.out.fcsgxresult
                        .map { meta, file ->
                            file ? [[id: meta.id], file] : [[:], []]
                        }

    emit:

    fcsgxresult
    genomedict         = samtools_reference
    fcsgx_report_txt
    fcsgx_taxonomy_rpt
    versions           = ch_versions

}
