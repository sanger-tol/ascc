// MODULE IMPORT BLOCK
include { BLAST_BLASTN                      } from '../../../modules/nf-core/blast/blastn/main'

include { SEQKIT_SLIDING                    } from '../../../modules/nf-core/seqkit/sliding/main'
include { BLAST_CHUNK_TO_FULL               } from '../../../modules/local/blast/chunk_to_full/main'
include { REFORMAT_FULL_OUTFMT6             } from '../../../modules/local/reformat/full_outfmt6/main'
include { BLAST_GET_TOP_HITS                } from '../../../modules/local/blast/get_top_hits/main'
include { GET_LINEAGE_FOR_TOP               } from '../../../modules/local/get/lineage_for_top/main'

workflow EXTRACT_NT_BLAST {
    take:
    input_genome            // Channel.of([ [ id: sample_id ], fasta ])
    blastn_db_path          // Channel.fromPath( db )
    ncbi_lineage_path       // Channel.fromPath( lineage_path )

    main:
    ch_versions             = Channel.empty()

    blastn_db_path
        .map { it ->
            [
                [id: "db"],
                it
            ]
        }
        .set { ch_blast }

    //
    // MODULE: CREATES A FASTA CONTAINING SLIDING WINDOWS OF THE INPUT GENOME
    //
    SEQKIT_SLIDING ( input_genome )
    ch_versions             = ch_versions.mix(SEQKIT_SLIDING.out.versions)

    //
    // MODULE: BLASTS THE INPUT GENOME AGAINST A LOCAL NCBI DATABASE
    //
    BLAST_BLASTN (
        SEQKIT_SLIDING.out.fastx,
        ch_blast
    )
    ch_versions             = ch_versions.mix(BLAST_BLASTN.out.versions)

    input_genome
        .map{ meta, file ->
            meta.id
        }
        .set { id }

    //
    // LOGIC: COLLECT THE BLAST OUTPUTS AND COLLECT THEM INTO ONE FILE
    //
    BLAST_BLASTN.out.txt
        .map { meta, files ->
            files
        }
        .collectFile(
            name: 'FULL_blast_results.txt',
            newLine: false
        )
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
    // MODULE: RETRIEVE LINEAGE INFORMATION FOR TOP BLAST RESULTS USING TAXIDS
    //
    GET_LINEAGE_FOR_TOP (
        BLAST_GET_TOP_HITS.out.tophits,
        ncbi_lineage_path
    )
    ch_versions             = ch_versions.mix(GET_LINEAGE_FOR_TOP.out.versions)

    // No conversion needed - BLAST results are already in the format expected by BlobToolKit

    emit:
    ch_blast_results        = BLAST_BLASTN.out.txt
    ch_formatted_results    = REFORMAT_FULL_OUTFMT6.out.full
    ch_top_lineages         = GET_LINEAGE_FOR_TOP.out.full
    ch_blast_hits           = BLAST_CHUNK_TO_FULL.out.full
    ch_btk_format           = BLAST_CHUNK_TO_FULL.out.full  // Format for BTK - full coordinates file
    versions                = ch_versions

}
