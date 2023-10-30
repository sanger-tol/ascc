#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { GET_KMER_COUNTS } from '../../modules/local/get_kmer_counts'
include { KMER_COUNT_DIM_REDUCTION } from '../../modules/local/kmer_count_dim_reduction'

workflow GET_KMERS_PROFILE {
    take:
    assembly_fasta                      // Channel [ val(meta), path(file) ]
    kmer_size                           // Channel [ val(integer) ]
    dimensionality_reduction_methods    // Channel [ val(string) ]
    n_neighbors_setting                 // Channel [ val(string) ]
    autoencoder_epochs_count            // Channel [ val(integer) ]

    main:
    ch_versions     = Channel.empty()

    assembly_fasta
        .map{ it ->
            tuple([id: it[0].id,
                    single_end: true
                ],
                it[1]
            )
        }
        .set { modified_input }

    //
    // MODULE: PRODUCE KMER COUNTS (USING KCOUNTER)
    //
    GET_KMER_COUNTS (
        modified_input,      // val(meta), path(reads)
        kmer_size            // val kmer_size
    )
    ch_versions = ch_versions.mix(GET_KMER_COUNTS.out.versions)


    //
    // LOGIC: CREATE CHANNEL OF LIST OF SELECTED METHODS
    //
    ch_methods = dimensionality_reduction_methods.splitCsv(sep: ',')





    //
    // MODULE: DIMENSIONALITY REDUCTION OF KMER COUNTS, USING SPECIFIED METHODS 
    //
    KMER_COUNT_DIM_REDUCTION (
        GET_KMER_COUNTS.out.csv,             // val(meta), path(kmer_counts)
        ch_methods,                         // val dimensionality_reduction_method
        n_neighbors_setting,                 // val n_neighbors_setting
        autoencoder_epochs_count            // val autoencoder_epochs_count
    )
    ch_versions = ch_versions.mix(KMER_COUNT_DIM_REDUCTION.out.versions)

    //
    // LOGIC: COMBINE OUTPUTS OF MULTIPLE METHODS
    //
    // ch_data
    //     .combine( alignment_datadir )
    //     .combine( assembly_classT )
    //     .map {
    //         ch_org, data_dir, classT ->
    //             file("${data_dir}${classT}/csv_data/${ch_org}-data.csv")
    //     }
    //     .splitCsv( header: true, sep:',')
    //     .map( row ->
    //     tuple([ org:    row.org,
    //             type:   row.type,
    //             id:     row.data_file.split('/')[-1].split('.MOD.')[0]
    //         ],
    //         file(row.data_file)
    //     ))
    //     .branch {
    //         pep: it[0].type  == 'pep'
    //         gen: it[0].type  == 'cdna'
    //         rna: it[0].type  == 'rna'
    //         cds: it[0].type  == 'cds'
    //     }
    //     .set {ch_alignment_data}



    emit:
    KMER_COUNT_DIM_REDUCTION.out.csv
    versions        = ch_versions.ifEmpty(null)
}
