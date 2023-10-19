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

    //
    // MODULE: RUN BLAST WITH GENOME AGAINST ORGANELLAR GENOME
    //
    BLAST_BLASTN (
        SED_SED.out.sed,
        BLAST_MAKEBLASTDB.out.db
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
        .branch {
            valid:      it[0].sz >= 250
            invalid:    it[0].sz <= 249
        }
        .set { blast_check }

    // How to get this printing blast_check.invalid.map{it[0].og} too?
    if ( blast_check.invalid ) {
        log.warn("BLAST RESULTS NOT AVAILABLE FOR $organellar_var.val")
    } else {
        log.warn("BLAST RESULTS FOUND FOR $organellar_var.val")
    }

    //
    // MODULE: FILTER COMMENTS OUT OF THE BLAST OUTPUT
    //
    FILTER_COMMENTS (
        blast_check.valid
    )
    ch_versions     = ch_versions.mix(FILTER_COMMENTS.out.versions)

    //
    // MODULE: EXTRACT CONTAMINANTS FROM THE BUSCO REPORT
    //
    EXTRACT_CONTAMINANTS (
        FILTER_COMMENTS.out.txt
    )
    ch_versions     = ch_versions.mix(EXTRACT_CONTAMINANTS.out.versions)

    //
    // LOGIC: COMBINE CHANNELS INTO FORMAT OF ID, ORGANELLE ID AND FILES
    //
    EXTRACT_CONTAMINANTS.out.bed
        .combine ( organellar_tuple )
        .map { blast_meta, blast_txt, organelle_meta, organelle_fasta ->
            tuple( [    id          :   blast_meta.id,
                        organelle   :   organelle_meta.id   ],
                    blast_txt
            )
        }
        .set { reformatted_recomendations }

    //
    // MODULE: GENERATE BED FILE OF ORGANELLAR SITES RECOMENDED TO BE REMOVED
    //
    ORGANELLE_CONTAMINATION_RECOMMENDATIONS (
        reformatted_recomendations
    )
    ch_versions     = ch_versions.mix(ORGANELLE_CONTAMINATION_RECOMMENDATIONS.out.versions)

    emit:
    organlle_report = ORGANELLE_CONTAMINATION_RECOMMENDATIONS.out.bed
    versions        = ch_versions.ifEmpty(null)

}
