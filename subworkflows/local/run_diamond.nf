include {   SEQKIT_SLIDING                                      } from '../../modules/nf-core/seqkit/sliding/main'
include {   DIAMOND_BLASTX                                      } from '../../modules/nf-core/diamond/blastx/main'
include {   BLAST_CHUNK_TO_FULL as DIAMOND_BLAST_CHUNK_TO_FULL  } from '../../modules/local/blast_chunk_to_full'
include {   REFORMAT_DIAMOND_OUTFMT6                            } from '../../modules/local/format_diamond_outfmt6'

workflow RUN_DIAMOND {
    take:
    reference_tuple     // tuple [[meta.id, meta.sliding, meta.window], reference]
    diamond_db          // val (path)

    main:
    ch_versions     = Channel.empty()
    ch_ext          = Channel.of(6)
    ch_columns      = Channel.of("qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms sphylums salltitles")

    //
    // MODULE: CREATE SLIDING WINDOW OF THE INPUT ASSEMBLY
    //
    SEQKIT_SLIDING (
        reference_tuple
    )
    ch_versions     = ch_versions.mix(SEQKIT_SLIDING.out.versions)

    //
    // MODULE: BLAST THE SLIDING WINDOW FASTA AGAINST THE DIAMOND DB
    //
    DIAMOND_BLASTX (
        SEQKIT_SLIDING.out.fastx,
        diamond_db,
        ch_ext,
        ch_columns
    )
    ch_versions     = ch_versions.mix(DIAMOND_BLASTX.out.versions)

    //
    // MODULE: COMBINE THE CHUNKS INTO THE FULL GENOME
    //
    DIAMOND_BLAST_CHUNK_TO_FULL (
        DIAMOND_BLASTX.out.txt
    )
    ch_versions     = ch_versions.mix(DIAMOND_BLAST_CHUNK_TO_FULL.out.versions)

    //
    // MODULE: REFORMAT THE DIAMOND OUTPUT
    //
    REFORMAT_DIAMOND_OUTFMT6 (
        DIAMOND_BLAST_CHUNK_TO_FULL.out.full
    )
    ch_versions     = ch_versions.mix(REFORMAT_DIAMOND_OUTFMT6.out.versions)

    emit:
    full            = DIAMOND_BLAST_CHUNK_TO_FULL.out.full
    reformed        = REFORMAT_DIAMOND_OUTFMT6.out.full
    versions        = ch_versions.ifEmpty(null)
}