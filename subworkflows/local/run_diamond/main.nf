//
// NF-CORE MODULE IMPORT
//
include { SEQKIT_SLIDING                                      } from '../../../modules/nf-core/seqkit/sliding/main'
include { DIAMOND_BLASTX                                      } from '../../../modules/nf-core/diamond/blastx/main'

//
// LOCAL MODULE IMPORTS
//
include { BLAST_CHUNK_TO_FULL as DIAMOND_BLAST_CHUNK_TO_FULL  } from '../../../modules/local/blast/chunk_to_full/main'
include { REFORMAT_TO_HITS_FILE                               } from '../../../modules/local/reformat/to_hits_file/main'
include { REFORMAT_DIAMOND_OUTFMT6                            } from '../../../modules/local/reformat/diamond_outfmt6/main'

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
    // LOGIC: GENERATE THE INPUT CHANNELS NEEDED FOR THE INPUT OF BLAST.
    //
    diamond_db
        .map{ it ->
            tuple(
                [id: "db"],
                it
            )
        }
        .set{diamond_db_path}

    SEQKIT_SLIDING.out.fastx
        .combine(ch_ext)
        .combine(ch_columns)
        .combine(diamond_db_path)
        .multiMap{ meta, reference, extensions, columns, meta2, db_path ->
            reference: tuple(meta, reference)
            db_path: tuple(meta2, db_path)
            ext_ch: extensions
            col_ch: columns
        }
        .set {blast_input}


    //
    // MODULE: BLAST THE SLIDING WINDOW FASTA AGAINST THE DIAMOND DB.
    //
    DIAMOND_BLASTX (
        blast_input.reference,
        blast_input.db_path,
        blast_input.ext_ch,
        blast_input.col_ch
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
    REFORMAT_TO_HITS_FILE(
        DIAMOND_BLAST_CHUNK_TO_FULL.out.full
    )
    ch_versions     = ch_versions.mix(REFORMAT_TO_HITS_FILE.out.versions)


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
    hits_file       = REFORMAT_TO_HITS_FILE.out.hits_file
    versions        = ch_versions
}
