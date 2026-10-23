include { SAMTOOLS_DICT         } from '../../../modules/nf-core/samtools/dict/main'
include { FCSGX_RUNGX           } from '../../../modules/sanger-tol/fcsgx/rungx/main'
include { FCSGX_PARSERESULTS    } from '../../../modules/sanger-tol/fcsgx/parseresults/main'

workflow FCSGX_PARSECSV {

    take:
    reference               // channel [ val(meta), path(file) ]
    fcsgxpath               // channel path(file)
    ncbi_rankedlineage_path // channel path(file)

    main:

    //
    // MODULE: USE SAMTOOLS_DICT TO GET THE ORIGIN FILE OF EACH SEQUENCE
    //         ITS NOT NEEDED BY ANYTHING IN THE SUBWORKFLOW BUT IS A NICE
    //         VERIFICATION OF SEQ ORIGIN
    //
    SAMTOOLS_DICT(
        reference
    )


    //
    // MODULE: RUN FCSGX FOR CLASSIFICATION OF SEQUENCES IN ASSEMBLY
    //
    FCSGX_RUNGX (
        reference.map { meta, ref -> tuple(meta, meta.taxid, ref) },
        fcsgxpath,
        [],
        "production" in workflow.profile.tokenize(',')
    )

    fcsgx_report_txt    = FCSGX_RUNGX.out.fcsgx_report
                            .map { meta, file ->
                                file ? tuple([ id: meta.id ], file) : [[:], []]
                            }

    fcsgx_taxonomy_rpt  = FCSGX_RUNGX.out.taxonomy_report
                            .map { meta, file ->
                                file ? tuple([ id: meta.id ], file) : [[:], []]
                            }


    //
    // MODULE: CONVER FCSGX_RUNGX RESULTS INTO A SINGLE CSV FILE
    //
    FCSGX_PARSERESULTS (
        fcsgx_taxonomy_rpt,
        fcsgx_report_txt,
        ncbi_rankedlineage_path
    )

    fcsgxresult     = FCSGX_PARSERESULTS.out.fcsgxresult
                        .map { meta, file ->
                            file ? tuple(meta, file) : [[:], []]
                        }

    emit:
    fcsgxresult
    genomedict         = SAMTOOLS_DICT.out.dict
    fcsgx_report_txt
    fcsgx_taxonomy_rpt
}
