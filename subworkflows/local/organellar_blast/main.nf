include { SED_SED                                    }   from '../../../modules/local/sed/sed/main'
include { BLAST_MAKEBLASTDB                          }   from '../../../modules/nf-core/blast/makeblastdb'
include { BLAST_BLASTN                               }   from '../../../modules/nf-core/blast/blastn'
include { EXTRACT_CONTAMINANTS                       }   from '../../../modules/local/extract/contaminants/main'
include { FILTER_COMMENTS                            }   from '../../../modules/local/filter/comments/main'
include { ORGANELLE_CONTAMINATION_RECOMMENDATIONS    }   from '../../../modules/local/organelle/contamination_recommendations/main'

//
// WORKFLOW: GENERATE A BED FILE CONTAINING LOCATIONS OF PUTATIVE ORGANELLAR SEQUENCE
//
workflow ORGANELLAR_BLAST {
    take:
    reference_tuple     // tuple([sample_id], reference_fasta)
    organellar_tuple    // tuple([organelle], organellar_fasta)

    main:
    ch_versions     = Channel.empty()

    reference_tuple
        .combine(organellar_tuple)
        .map { ref_meta, ref_file, org_meta, org_file ->
            def meta = ref_meta + [ og: org_meta.id]
            tuple(meta, ref_file)
        }
        .set{ new_ref_tuple }


    //
    // MODULE: STRIP SPACES OUT OF GENOMIC FASTA
    //
    SED_SED (
        new_ref_tuple
    )
    ch_versions     = ch_versions.mix(SED_SED.out.versions)


    //
    // MODULE: GENERATE BLAST DB ON ORGANELLAR GENOME
    //
    BLAST_MAKEBLASTDB (
        organellar_tuple
    )
    ch_versions     = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)


    //
    // MODULE: RUN BLAST WITH GENOME AGAINST ORGANELLAR GENOME
    //
    SED_SED.out.sed
        .combine(BLAST_MAKEBLASTDB.out.db)
        .multiMap{ meta, ref, meta2, blast_db ->
            reference_tuple: tuple(meta, ref)
            blastdb_tuple: tuple(meta, blast_db)
        }
        .set { ref_and_db }

    BLAST_BLASTN (
        ref_and_db.reference_tuple,
        ref_and_db.blastdb_tuple
    )
    ch_versions     = ch_versions.mix(BLAST_BLASTN.out.versions)


    //
    // LOGIC: REORGANISE CHANNEL FOR DOWNSTREAM PROCESS
    //
    BLAST_BLASTN.out.txt
        .combine ( organellar_tuple )
        .map { meta, file, org_meta, org_file ->
            tuple ( [   id: meta.id,                // Assembly Name
                        og: org_meta.id,            // Organellar Name
                        sz: file.size() ],          // Size of assembly
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

    no_comments.valid.view{"ORGANELLAR_BLAST VALID ORGANELLAR: $it"}

    // Strip out a ton of junk meta
    new_ref_tuple
        .map{meta, file ->
            [[id: meta.id, og: meta.og], file]
        }
        .set{fixed_ref}

    fixed_ref.view{"ORGANELLAR_BLAST REFERENCE: $it"}


    no_comments
        .valid
        .map{ meta, file ->
            [[id: meta.id, og: meta.og], file]
        }
        .combine(fixed_ref, by: 0)
        .multiMap { meta, no_comment_file, reference ->
            filtered: tuple(meta, no_comment_file)
            reference: tuple(meta, reference)

        }
        .set { mapped }

    mapped.filtered.view{"FILTERED: $it"}
    mapped.reference.view{"REF: $it"}


    //
    // MODULE: EXTRACT CONTAMINANTS FROM THE BLAST REPORT
    //
    EXTRACT_CONTAMINANTS (
        mapped.filtered,
        mapped.reference
    )
    ch_versions     = ch_versions.mix(EXTRACT_CONTAMINANTS.out.versions)


    //
    // LOGIC: COMBINE CHANNELS INTO FORMAT OF ID, ORGANELLE ID AND FILES
    //
    EXTRACT_CONTAMINANTS.out.contamination_bed
        .combine ( organellar_tuple)
        .map { blast_meta, blast_txt, organelle_meta, organelle_fasta ->
            tuple( [    id          :   blast_meta.id,
                        organelle   :   organelle_meta.id   ],
                    file(blast_txt)
            )
        }
        .set { reformatted_recommendations }


    //
    // MODULE: GENERATE BED FILE OF ORGANELLAR SITES RECOMENDED TO BE REMOVED
    //
    ORGANELLE_CONTAMINATION_RECOMMENDATIONS (
        reformatted_recommendations
    )
    ch_versions     = ch_versions.mix(ORGANELLE_CONTAMINATION_RECOMMENDATIONS.out.versions)


    emit:
    organelle_report= ORGANELLE_CONTAMINATION_RECOMMENDATIONS.out.recommendations
    versions        = ch_versions

}
