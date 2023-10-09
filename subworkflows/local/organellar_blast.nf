import { SED_SED                                    }   from '../modules/local/sed_sed'
import { BLAST_MAKEBLASTDB                          }   from '../modules/nf-core/blast/makeblastdb'
import { BLAST_BLASTN                               }   from '../modules/nf-core/blast/blastn'
import { EXTRACT_CONTAMINANTS                       }   from '../modules/local/extract_contaminants'
import { FILTER_COMMENTS                            }   from '../modules/local/filter_comments'
import { ORGANELLE_CONTAMINATION_RECOMMENDATIONS    }   from '../modules/local/organellar_contamination_report'

//
// WORKFLOW: GENERATE A BED FILE CONTAINING LOCATIONS OF PUTATIVE ORGANELLAR SEQUENCE
//
workflow ORGANELLAR_BLAST {
    take:
    reference_tuple     // tuple([sample_id], reference_fasta)
    organellar_tuple    // tuple([organelle], organellar_fasta)

    main:

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
    BLAST_MAKEBLASTDB (
        organellar_tuple
    )
    ch_versions     = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

    //
    // MODULE: RUN BLAST WITH GENOME AGAINST ORGANELLAR GENOME
    //
    BLAST_BLASTN (
        SED_SED.out.sed,
        BLAST_MAKEBLASTDB.out.db
    )
    ch_versions     = ch_versions.mix(BLASTN.out.versions)

    //
    // MODULE: FILTER COMMENTS OUT OF THE BLAST OUTPUT
    //
    FILTER_COMMENTS (
        BLAST_BLASTN.out.txt
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
    // MODULE: GENERATE BED FILE OF ORGANELLAR SITES RECOMENDED TO BE REMOVED
    //
    ORGANELLE_CONTAMINATION_RECOMMENDATIONS (
        EXTRACT_CONTAMINANTS.out.bed
    )
    ch_versions     = ch_versions.mix(ORGANELLA_CONTAM_REPORT.out.versions)

    emit:
    organlle_report = ORGANELLE_CONTAMINATION_RECOMMENDATIONS.out.bed
    versions        = ch_versions.ifEmpty(null)

}
