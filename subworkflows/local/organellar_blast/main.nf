//
// NF-CORE MODULE IMPORT
//
include { BLAST_MAKEBLASTDB                          }   from '../../../modules/nf-core/blast/makeblastdb'
include { BLAST_BLASTN                               }   from '../../../modules/nf-core/blast/blastn'

//
// LOCAL MODULE IMPORT
//
include { SED_SED                                    }   from '../../../modules/local/sed/sed/main'
include { EXTRACT_CONTAMINANTS                       }   from '../../../modules/local/extract/contaminants/main'
include { FILTER_COMMENTS                            }   from '../../../modules/local/filter/comments/main'
include { ORGANELLE_CONTAMINATION_RECOMMENDATIONS    }   from '../../../modules/local/organelle/contamination_recommendations/main'

//
// WORKFLOW: GENERATE A BED FILE CONTAINING LOCATIONS OF PUTATIVE ORGANELLAR SEQUENCE
//
workflow ORGANELLAR_BLAST {
    take:
    reference_tuple     // channel [sample_id], reference_fasta
    organellar_tuple    // channel [organelle], organellar_fasta

    main:
    ch_versions     = channel.empty()

    combined_refs = reference_tuple
        .combine(organellar_tuple)
        .multiMap { ref_meta, ref_file, org_meta, org_file ->
            def meta = ref_meta + [og: org_meta.id]
            refs: [meta, ref_file]
            orgs: [meta, org_file]
        }


    organellar_dbs_input = combined_refs.orgs
        .map { meta, org_file -> [[id: meta.og], org_file] }
        .unique { meta, file -> meta.id }

    //
    // MODULE: GENERATE BLAST DB ON ORGANELLAR GENOME
    //
    BLAST_MAKEBLASTDB (
        organellar_dbs_input,
        []
    )


    //
    // MODULE: RUN BLAST WITH GENOME AGAINST ORGANELLAR GENOME
    //
    ref_and_db = combined_refs.refs
        .map { meta, ref -> [meta.og, meta, ref] }
        .combine(
            BLAST_MAKEBLASTDB.out.db
                .map { meta, db -> [meta.id, db] },
            by: 0
        )
        .map { og, meta, ref, db -> [meta, ref, db] }
        .multiMap{ meta, ref, blast_db ->
            reference_tuple:    [meta, ref]
            blastdb_tuple:      [meta, blast_db]
        }


    BLAST_BLASTN (
        ref_and_db.reference_tuple,
        ref_and_db.blastdb_tuple,
        [],
        [],
        []
    )


    //
    // LOGIC: REORGANISE CHANNEL FOR DOWNSTREAM PROCESS
    //
    BLAST_BLASTN.out.txt
        .map { meta, file ->
            [[  id: meta.id,                // Assembly Name
                og: meta.og,                // Organellar Name (already in meta)
                sz: file.size()             // Size of assembly
            ], file ]
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
            def lines_in_file = file.countLines()
            valid: lines_in_file >= 1
            invalid : lines_in_file < 1

            log.info "[ASCC INFO] ORGANELLAR_BLAST results contain ${ lines_in_file } lines (> 0 is Valid)"
            log.info "\t--$meta.id & $meta.og "
        }
        .set { no_comments }


    //
    // NOTE: Strip out a ton of junk meta so we can join the tuples together
    //
    no_comments
        .valid
        .map{ meta, file ->
            tuple([id: meta.id, og: meta.og], file)
        }
        .combine(
            combined_refs.refs
                .map{ meta, file ->
                    tuple([id: meta.id, og: meta.og], file)
                },
            by: 0
        )
        .multiMap { meta, no_comment_file, assembly ->
            filtered: tuple(meta, no_comment_file)
            reference_ch: tuple(meta, assembly)

        }
        .set { mapped }


    //
    // MODULE: EXTRACT CONTAMINANTS FROM THE BLAST REPORT
    //
    EXTRACT_CONTAMINANTS (
        mapped.filtered,
        mapped.reference_ch
    )
    ch_versions     = ch_versions.mix(EXTRACT_CONTAMINANTS.out.versions)


    //
    // LOGIC: REFORMAT OUTPUT WITH ID AND ORGANELLE ID
    // NOTE: og is already in metadata from the initial combine flow.
    //
    reformatted_recommendations = EXTRACT_CONTAMINANTS.out.contamination_bed
        .map { meta, blast_txt ->
            tuple(
                [   id:         meta.id,
                    organelle:  meta.og
                ],
                blast_txt
            )
        }


    //
    // MODULE: GENERATE BED FILE OF ORGANELLAR SITES RECOMENDED TO BE REMOVED
    //
    ORGANELLE_CONTAMINATION_RECOMMENDATIONS (
        reformatted_recommendations
    )
    ch_versions     = ch_versions.mix(ORGANELLE_CONTAMINATION_RECOMMENDATIONS.out.versions)


    emit:
    organelle_report        = ORGANELLE_CONTAMINATION_RECOMMENDATIONS.out.recommendations
    full_organelle_report   = reformatted_recommendations
    versions                = ch_versions

}
