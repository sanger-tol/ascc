include {   SEQKIT_SLIDING                                      } from '../../modules/nf-core/seqkit/sliding/main'
include {   DIAMOND_BLASTX                                      } from '../../modules/nf-core/diamond/blastx/main'
include {   BLAST_CHUNK_TO_FULL as DIAMOND_BLAST_CHUNK_TO_FULL  } from '../../modules/local/blast_chunk_to_full'
include {   CONVERT_TO_HITS_FILE                                } from '../../modules/local/convert_to_hits_file'
include {   REFORMAT_DIAMOND_OUTFMT6                            } from '../../modules/local/format_diamond_outfmt6'

workflow RUN_DIAMOND {
    take:
    reference_tuple     // tuple [[meta.id, meta.sliding, meta.window], reference]
    diamond_db          // val (path)

    main:
    ch_versions     = Channel.empty()
    ch_ext          = Channel.of("txt")
    ch_columns      = Channel.of("qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms sphylums salltitles")

    //
    // MODULE: CREATE SLIDING WINDOW OF THE INPUT ASSEMBLY
    //
    SEQKIT_SLIDING (
        reference_tuple
    )
    ch_versions     = ch_versions.mix(SEQKIT_SLIDING.out.versions)

    //
    // LOGIC: 
    //
    SEQKIT_SLIDING.out.fastx
        .combine(diamond_db)
        .combine(ch_ext)
        .combine(ch_columns)
        .multiMap{meta_slide, files_slide, diamond, file_ext, blast_columns ->
            seqkit_out: tuple(meta_slide, files_slide)
            diamond_db: diamond
            extension: file_ext
            columns: blast_columns
        }
        .set{blast_input}

    //
    // MODULE: BLAST THE SLIDING WINDOW FASTA AGAINST THE DIAMOND DB
    //
    DIAMOND_BLASTX (
        blast_input.seqkit_out,
        blast_input.diamond_db,
        blast_input.extension,
        blast_input.columns
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
    // MODULE: CONVERT THE FULL GENOME FILE INTO A HITS FILE
    //
    CONVERT_TO_HITS_FILE(
        DIAMOND_BLAST_CHUNK_TO_FULL.out.full
    )
    ch_versions     = ch_versions.mix(CONVERT_TO_HITS_FILE.out.versions)

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
    hits_file       = CONVERT_TO_HITS_FILE.out.hits_file
    versions        = ch_versions.ifEmpty(null)
}