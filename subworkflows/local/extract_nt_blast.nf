// MODULE IMPORT BLOCK
include { BLAST_BLASTN          } from '../../modules/nf-core/blast/blastn/main'

include { SEQKIT_SLIDING        } from '../../modules/nf-core/seqkit/sliding/main'
include { BLAST_CHUNK_TO_FULL   } from '../../modules/local/blast_chunk_to_full'
include { REFORMAT_FULL_OUTFMT6 } from '../../modules/local/reformat_full_outfmt6'
include { BLAST_GET_TOP_HITS    } from '../../modules/local/blast_get_top_hits'
include { GET_LINEAGE_FOR_TOP   } from '../../modules/local/get_lineage_for_top'

workflow EXTRACT_NT_BLAST {
    take:
    input_genome            // Channel.of([ [ id: sample_id ], fasta ])
    blastn_db_path          // Channel.of( path )
    ncbi_taxonomy_path      // Channel.of( path )
    ncbi_lineage_path       // Channel.of( path )

    main:
    ch_versions             = Channel.empty()

    //
    // MODULE: CREATES A FASTA CONTAINING SLIDING WINDOWS OF THE INPUT GENOME
    // NOTE: needs to be installed once accepted by NF-core
    SEQKIT_SLIDING ( input_genome )
    ch_versions             = ch_versions.mix(SEQKIT_SLIDING.out.versions)

    //
    // MODULE: BLASTS THE INPUT GENOME AGAINST A LOCAL NCBI DATABASE
    //
    BLAST_BLASTN (
        SEQKIT_SLIDING.out.fastx,
        blastn_db_path
    )
    ch_versions             = ch_versions.mix(BLAST_BLASTN.out.versions)

    //
    // LOGIC: MERGE WITH INPUT GENOME TO TAKE SAMPLE_ID
    //
    BLAST_BLASTN.out.txt
        .combine ( input_genome )
        .map { it ->
            tuple( it[1],
                    it[2],
                    it[0]
            )
        }
        .set { mapped_db }

    //
    // MODULE:
    //
    BLAST_CHUNK_TO_FULL ( mapped_db )
    ch_versions             = ch_versions.mix(BLAST_CHUNK_TO_FULL.out.versions)

    //
    // MODULE:
    //
    REFORMAT_FULL_OUTFMT6 ( BLAST_CHUNK_TO_FULL.out.full )
    ch_versions             = ch_versions.mix(REFORMAT_FULL_OUTFMT6.out.versions)

    //
    // MODULE:
    //
    BLAST_GET_TOP_HITS ( REFORMAT_FULL_OUTFMT6.out.full )
    ch_versions             = ch_versions.mix(BLAST_GET_TOP_HITS.out.versions)

    //
    // MODULE:
    //
    GET_LINEAGE_FOR_TOP (
        BLAST_GET_TOP_HITS.out.tophits,
        ncbi_taxonomy_path,
        ncbi_lineage_path
    )
    ch_versions             = ch_versions.mix(GET_LINEAGE_FOR_TOP.out.versions)


    emit:
    versions                = ch_versions.ifEmpty(null)

}
