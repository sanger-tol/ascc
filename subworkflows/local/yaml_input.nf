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
                pacbio_reads:                                   ( data.pacbio_reads_path                )
                assembly_path:                                  ( file(data.assembly_path)              )
                pacbio_multiplexing_barcode_names:              ( data.pacbio_multiplexing_barcode_names)
                sci_name:                                       ( data.sci_name                         )
                taxid:                                          ( data.taxid                            )
                mito_fasta_path:                                ( data.mito_fasta_path                  )
                plastid_fasta_path:                             ( data.plastid_fasta_path               )
                nt_database:                                    ( data.nt_database                      )
                reference_proteomes:                            ( data.reference_proteomes              )
                nt_kraken_db_path:                              ( data.nt_kraken_db_path                )
                kmer_len:                                       ( data.kmer_len                         )
                fcs_gx_database_path:                           ( data.fcs_gx_database_path             )
                ncbi_taxonomy_path:                             ( data.ncbi_taxonomy_path               )
                ncbi_rankedlineage_path:                        ( data.ncbi_rankedlineage_path          )
                busco_lineages_folder:                          ( data.busco_lineages_folder            )
                seqkit_values:                                  ( data.seqkit                           )
                diamond_uniprot_database_path:                  ( data.diamond_uniprot_database_path    )
                diamond_nr_database_path:                       ( data.diamond_nr_database_path         )
                vecscreen_database_path:                        ( data.vecscreen_database_path          )

        }
        .set{ group }

    group
        .seqkit_values
        .flatten()
        .multiMap { data ->
            sliding_value                           :           ( data.sliding                          )
            window_value                            :           ( data.window                           )
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
        .combine( group.pacbio_reads )
        .map { id, file ->
            tuple(  [   id: id ],
                    file
            )
        }
        .set { ch_pacbio }

    Channel
        .fromPath( group.nt_database, checkIfExists=true )
        .set { blast_db }

    emit:
    reference_tuple                  = ch_reference
    pacbio_tuple                     = ch_pacbio
    assembly_title                   = group.assembly_title
    taxid                            = group.taxid
    nt_database                      = blast_nt_db
    nt_kraken_db_path                = group.nt_kraken_db_path
    ncbi_taxonomy_path               = group.ncbi_taxonomy_path
    ncbi_rankedlineage_path          = group.ncbi_rankedlineage_path
    busco_lineages_folder            = group.busco_lineages_folder
    fcs_gx_database_path             = group.fcs_gx_database_path
    diamond_uniprot_database_path    = group.diamond_uniprot_database_path
    diamond_nr_database_path         = group.diamond_nr_database_path
    vecscreen_database_path          = group.vecscreen_database_path
    seqkit_sliding                   = seqkit.sliding_value
    seqkit_window                    = seqkit.window_value
    versions                         = ch_versions.ifEmpty(null)
}

def readYAML( yamlfile ) {
    return new Yaml().load( new FileReader( yamlfile.toString() ) )
}
