#!/usr/bin/env nextflow

import org.yaml.snakeyaml.Yaml

workflow YAML_INPUT {
    take:
    input_file  // input_yaml_from_commandline

    main:
    ch_versions = Channel.empty()

    input_file
        .map { file -> readYAML(file) }
        .set { yamlfile }

    //
    // LOGIC: PARSES THE TOP LEVEL OF YAML VALUES
    //
    yamlfile
        .flatten()
        .multiMap { data ->
                assembly_title:                                 ( data.assembly_title                   )
                reads_path:                                     ( data.reads_path                       )
                reads_type:                                     ( data.reads_type                       )
                assembly_path:                                  ( file(data.assembly_path)              )
                pacbio_barcodes:                                ( file(data.pacbio_barcodes)            )
                pacbio_multiplexing_barcode_names:              ( data.pacbio_multiplexing_barcode_names)
                sci_name:                                       ( data.sci_name                         )
                taxid:                                          ( data.taxid                            )
                mito_fasta_path:                                ( data.mito_fasta_path                  )
                plastid_fasta_path:                             ( data.plastid_fasta_path               )
                nt_database:                                    ( data.nt_database                      )
                reference_proteomes:                            ( data.reference_proteomes              )
                nt_kraken_db_path:                              ( data.nt_kraken_db_path                )
                kmer_len:                                       ( data.kmer_len                         )
                dimensionality_reduction_methods:               ( data.dimensionality_reduction_methods )
                fcs_gx_database_path:                           ( data.fcs_gx_database_path             )
                ncbi_taxonomy_path:                             ( data.ncbi_taxonomy_path               )
                ncbi_rankedlineage_path:                        ( data.ncbi_rankedlineage_path          )
                ncbi_accessionids:                              ( data.ncbi_accessionids_folder         )
                busco_lineages_folder:                          ( data.busco_lineages_folder            )
                seqkit_values:                                  ( data.seqkit                           )
                diamond_uniprot_database_path:                  ( data.diamond_uniprot_database_path    )
                diamond_nr_database_path:                       ( data.diamond_nr_database_path         )
                vecscreen_database_path:                        ( data.vecscreen_database_path          )
                neighbours:                                     ( data.n_neighbours                     )
        }
        .set{ group }

    group
        .seqkit_values
        .flatten()
        .multiMap { data ->
            sliding_value                           :           ( data.sliding                                  )
            window_value                            :           ( data.window                                   )
        }
        .set { seqkit }

    group.assembly_title
        .combine( group.assembly_path )
        .map { id, file ->
            tuple(  [   id: id ],
                    file
            )
        }
        .set { ch_reference }

    group.assembly_title
        .combine( group.reads_path )
        .map { id, file ->
            tuple(  [   id: id ],
                    file
            )
        }
        .set { ch_pacbio }

    group.pacbio_barcodes
        .map { file ->
            tuple(  [   id: "pacbio barcodes"   ],
                    file
            )
        }
        .set { ch_barcodes }

    group.mito_fasta_path
        .map{ it ->
            tuple(  [   id: "mitochondrial_genome"  ],
                    it
            )
        }
        .set{ ch_mito }

    group.plastid_fasta_path
        .map{ it ->
            tuple(  [   id: "plastid_genome"    ],
                    it
            )
        }
        .set{ ch_plastid }

    emit:
    reference_tuple                  = ch_reference
    pacbio_tuple                     = ch_pacbio
    reads_type                       = group.reads_type
    pacbio_barcodes                  = ch_barcodes
    pacbio_multiplex_codes           = group.pacbio_multiplexing_barcode_names
    assembly_title                   = group.assembly_title
    assembly_path                    = group.assembly_path
    taxid                            = group.taxid
    nt_database                      = group.nt_database
    nt_kraken_db_path                = group.nt_kraken_db_path
    ncbi_accessions                  = group.ncbi_accessionids
    ncbi_taxonomy_path               = group.ncbi_taxonomy_path
    ncbi_rankedlineage_path          = group.ncbi_rankedlineage_path
    busco_lineages_folder            = group.busco_lineages_folder
    fcs_gx_database_path             = group.fcs_gx_database_path
    diamond_uniprot_database_path    = group.diamond_uniprot_database_path
    diamond_nr_database_path         = group.diamond_nr_database_path
    vecscreen_database_path          = group.vecscreen_database_path
    seqkit_sliding                   = seqkit.sliding_value
    seqkit_window                    = seqkit.window_value
    dimensionality_reduction_methods = group.dimensionality_reduction_methods
    mito_tuple                       = ch_mito
    mito_var                         = "mitochondrial_genome"
    plastid_tuple                    = ch_plastid
    plastid_var                      = "plastid_genome"
    kmer_len                         = group.kmer_len
    n_neighbours                     = group.neighbours
    versions                         = ch_versions.ifEmpty(null)
}

def readYAML( yamlfile ) {
    return new Yaml().load( new FileReader( yamlfile.toString() ) )
}
