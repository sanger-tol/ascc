// MODULE IMPORT BLOCK
include { BLAST_V5_DATABASE     } from '../../modules/local/blast_v5_database'
include { BLAST_BLASTN as BLAST_BLASTN_MOD }   from '../../modules/nf-core/blast/blastn'

include { SEQKIT_SLIDING        } from '../../modules/nf-core/seqkit/sliding/main'
include { BLAST_CHUNK_TO_FULL   } from '../../modules/local/blast_chunk_to_full'
include { REFORMAT_FULL_OUTFMT6 } from '../../modules/local/reformat_full_outfmt6'
include { BLAST_GET_TOP_HITS    } from '../../modules/local/blast_get_top_hits'
include { GET_LINEAGE_FOR_TOP   } from '../../modules/local/get_lineage_for_top'

workflow EXTRACT_NT_BLAST {
    take:
    input_genome            // Channel.of([ [ id: sample_id ], fasta ])
    blastn_db_path          // Channel.of( path )
    ncbi_accessions         // Channel.of( path )
    ncbi_lineage_path       // Channel.of( path )

    main:
    ch_versions             = Channel.empty()

    //
    // MODULE: CREATES A FASTA CONTAINING SLIDING WINDOWS OF THE INPUT GENOME
    //
    SEQKIT_SLIDING ( input_genome )
    ch_versions             = ch_versions.mix(SEQKIT_SLIDING.out.versions)

    //
    // MODULE: BLASTS THE INPUT GENOME AGAINST A LOCAL NCBI DATABASE
    //
    BLAST_BLASTN_MOD (
        SEQKIT_SLIDING.out.fastx,
        [[id: "db"], blastn_db_path]
    )
    ch_versions             = ch_versions.mix(BLAST_BLASTN_MOD.out.versions)

    input_genome
        .map{ meta, file ->
            meta.id
        }
        .set { id }

    //
    // LOGIC: COLLECT THE BLAST OUTPUTS AND COLLECT THEM INTO ONE FILE
    //
    BLAST_BLASTN_MOD.out.txt
        .map { meta, files ->
            files
        }
        .collectFile( name: 'FULL_blast_results.txt', newLine: false)
        .combine( id )
        .map { file, identity ->
            tuple(  [   id: identity    ],
                    file
                )
            }
        .set { blast_results }

    //
    // MODULE: CONVERT CHUNK_COORDINATES TO FULL_COORINDATES
    //
    BLAST_CHUNK_TO_FULL ( blast_results )
    ch_versions             = ch_versions.mix(BLAST_CHUNK_TO_FULL.out.versions)


    //
    // MODULE: RE_ORDER THE DATA IN THE FULL_COORDINATE FILE
    //
    REFORMAT_FULL_OUTFMT6 ( BLAST_CHUNK_TO_FULL.out.full )
    ch_versions             = ch_versions.mix(REFORMAT_FULL_OUTFMT6.out.versions)


    //
    // LOGIC: BRANCH DEPENDING ON WHETHER FILE HAS CONTENTS
    //
    REFORMAT_FULL_OUTFMT6.out.full
        .branch {
            valid:      it[1].readLines().size() >= 1
            invalid:    true
        }
        .set { gatekeeper }


    //
    // MODULE: RETURN ONLY THE TOP HITS PER SEQUENCE
    //
    BLAST_GET_TOP_HITS (
        gatekeeper.valid
    )
    ch_versions             = ch_versions.mix(BLAST_GET_TOP_HITS.out.versions)


    //
    // MODULE: USING THE accession2taxid DATABASE WE RETREIVE LINEAGE INFORMATION PER TOP RESULT
    //
    GET_LINEAGE_FOR_TOP (
        BLAST_GET_TOP_HITS.out.tophits,
        ncbi_accessions,
        ncbi_lineage_path
    )
    ch_versions             = ch_versions.mix(GET_LINEAGE_FOR_TOP.out.versions)


    emit:
    ch_top_lineages         = GET_LINEAGE_FOR_TOP.out.full
    ch_blast_hits           = BLAST_CHUNK_TO_FULL.out.full
    versions                = ch_versions.ifEmpty(null)

}

process get_string {
    input:
    val(nin)

    output:
    stdout

    script:
    """
    echo $nin
    """
}
