include { SED_SED                                    }   from '../../modules/local/sed_sed'
include { BLAST_MAKEBLASTDB                          }   from '../../modules/nf-core/blast/makeblastdb'
include { BLAST_BLASTN                               }   from '../../modules/nf-core/blast/blastn'
include { EXTRACT_CONTAMINANTS                       }   from '../../modules/local/extract_contaminants'
include { FILTER_COMMENTS                            }   from '../../modules/local/filter_comments'
include { ORGANELLE_CONTAMINATION_RECOMMENDATIONS    }   from '../../modules/local/organelle_contamination_recommendations'

//
// WORKFLOW: GENERATE A BED FILE CONTAINING LOCATIONS OF PUTATIVE ORGANELLAR SEQUENCE
//
workflow ORGANELLAR_BLAST {
    take:
    reference_tuple     // tuple([sample_id], reference_fasta)
    organellar_var      // str
    organellar_tuple    // tuple([organelle], organellar_fasta)

    main:
    ch_versions             = Channel.empty()

    //
    // MODULE: STRIP SPACES OUT OF GENOMIC FASTA
    //
    SED_SED (
        reference_tuple
    )
    ch_versions     = ch_versions.mix(SED_SED.out.versions)

    //
    // MODULE: GENERATE BLAST DB ON ORGANELLAR GENOME
    //
    organellar_tuple
        .map{it ->
            it[1]
        }
        .set { organelle_file }

    BLAST_MAKEBLASTDB (
        organelle_file
    )
    ch_versions     = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

    BLAST_MAKEBLASTDB.out.db
        .combine( organellar_tuple )
        .map { db_path, meta, organelle ->
            tuple ( meta,
                    db_path
            )
        }
        .set { organellar_check_db }

    //
    // MODULE: RUN BLAST WITH GENOME AGAINST ORGANELLAR GENOME
    //
    BLAST_BLASTN (
        SED_SED.out.sed,
        organellar_check_db
    )
    ch_versions     = ch_versions.mix(BLAST_BLASTN.out.versions)

    //
    // LOGIC: FILTER BLAST RESULTS WITH ONLY COMMENT LINES THESE SHOULD BE THOSE UNDER 250bytes
    //
    BLAST_BLASTN.out.txt
        .combine ( organellar_tuple )
        .map { meta, file, org_meta, org_file ->
            tuple ( [   id: meta.id,
                        og: org_meta.id,
                        sz: file.size() ],
                    file
            )
        }
        .set { blast_check }

    //
    // MODULE: FILTER COMMENTS OUT OF THE BLAST OUTPUT, ALSO BLAST result
    //
    FILTER_COMMENTS (
        blast_check
    )
    ch_versions     = ch_versions.mix(FILTER_COMMENTS.out.versions)

    //
    // LOGIC: IF FILTER_COMMENTS RETURNS FILE WITH NO LINES THEN SUBWORKFLOWS STOPS
    //
    FILTER_COMMENTS.out.txt
        .branch { meta, file ->
            valid: file.countLines() >= 1
            invalid : file.countLines() < 1
        }
        .set {no_comments}

    //
    // MODULE: EXTRACT CONTAMINANTS FROM THE BLAST REPORT
    //
    EXTRACT_CONTAMINANTS (
        FILTER_COMMENTS.out.txt,
        reference_tuple
    )
    ch_versions     = ch_versions.mix(EXTRACT_CONTAMINANTS.out.versions)

    //
    // LOGIC: COMBINE CHANNELS INTO FORMAT OF ID, ORGANELLE ID AND FILES
    //
    EXTRACT_CONTAMINANTS.out.contamination_bed
        .combine ( organellar_tuple )
        .map { blast_meta, blast_txt, organelle_meta, organelle_fasta ->
            tuple( [    id          :   blast_meta.id,
                        organelle   :   organelle_meta.id   ],
                    blast_txt
            )
        }
        .set { reformatted_recomendations }

    reformatted_recomendations.view()

    //
    // MODULE: GENERATE BED FILE OF ORGANELLAR SITES RECOMENDED TO BE REMOVED
    //
    ORGANELLE_CONTAMINATION_RECOMMENDATIONS (
        reformatted_recomendations
    )
    ch_versions     = ch_versions.mix(ORGANELLE_CONTAMINATION_RECOMMENDATIONS.out.versions)

    emit:
    organelle_report = ORGANELLE_CONTAMINATION_RECOMMENDATIONS.out.recomendations
    versions        = ch_versions.ifEmpty(null)

}
