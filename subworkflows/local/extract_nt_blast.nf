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
    // LOGIC: GLOB ALL *NIN FILES IN DIRECTORY AND SPLIT INTO CHANNELS
    //

    blastn_db_path
        .map(
            it -> file("${it}*") // glob all files in directory
        )
        .flatten() // flatten to file per channel
        .map(
            it ->
                tuple (
                    [   id: it.toString().split('/')[-1].split("\\....\$")[0] ], // get basename and trim off the extension, returns database prefix
                    it                                                           // list of files
                )
        )
        .groupTuple() // group files by id (which = db prefix)
        .map {
            meta, files ->
                tuple (
                    [   id: meta.id,
                        file_count: files.size() ], // get number of files
                        files
                )
        }
        .filter { it[0].file_count >= 8 } // a database is made of 8 files, less than this means it is an accessory to the db
        .set { databases_by_prefix }

    databases_by_prefix
        .combine( blastn_db_path )
        .map { meta, files, rootpath ->
            tuple( rootpath, meta.id )
        }
        .combine ( SEQKIT_SLIDING.out.fastx )
        .multiMap { root, db_prefix, meta, ref ->
            reference:  tuple( [ id:     meta.id ],
                                ref
                            )
            nin_db:     tuple( [ id:    db_prefix   ],
                                root
                            )
        }
        .set { nin }

    //
    // MODULE: BLASTS THE INPUT GENOME AGAINST A LOCAL NCBI DATABASE
    //
    BLAST_BLASTN_MOD (
        nin.reference,
        nin.nin_db
    )
    ch_versions             = ch_versions.mix(BLAST_BLASTN_MOD.out.versions)

    input_genome
        .map{ meta, file ->
            meta.id
        }
        .set { id }

    BLAST_BLASTN_MOD.out.txt
        .map { meta, files ->
            files
        }
        .collectFile( name: 'FULL_blast_results.txt', newLine: false) // concats all input files into one file!
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
        .map { meta, file ->
            tuple(  [   id: meta.id,
                        sz: file.size() ],
                    file
            )
        }
        .branch {
            valid:      it[0].sz >= 1
            invalid:    it[0].sz <= 0
        }
        .set { gatekeeper }

    //
    // MODULE:
    //
    BLAST_GET_TOP_HITS (
        gatekeeper.valid
    )
    ch_versions             = ch_versions.mix(BLAST_GET_TOP_HITS.out.versions)

    //
    // MODULE:
    //
    GET_LINEAGE_FOR_TOP (
        BLAST_GET_TOP_HITS.out.tophits,
        ncbi_accessions,
        ncbi_lineage_path
    )
    ch_versions             = ch_versions.mix(GET_LINEAGE_FOR_TOP.out.versions)

    emit:
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
