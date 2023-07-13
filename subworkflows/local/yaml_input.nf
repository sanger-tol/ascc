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
                assembly_title:                                 ( data.assembly_title )
                pacbio_reads:                                   ( data.pacbio_reads_path )
                assembly_path:                                  ( file(data.assembly_path) )
                pacbio_multiplexing_barcode_names:              ( data.pacbio_multiplexing_barcode_names )
                sci_name:                                       ( data.sci_name )
                taxid:                                          ( data.taxid )
                mito_fasta_path:                                ( data.mito_fasta_path )
                plastid_fasta_path:                             ( data.plastid_fasta_path )
                nt_database:                                    ( data.nt_database )
                reference_proteomes:                            ( data.reference_proteomes )
                nt_kraken_db_path:                              ( data.nt_kraken_db_path )
                kmer_len:                                       ( data.kmer_len )
                ncbi_taxonomy_path:                             ( data.ncbi_taxonomy_path )
                ncbi_lineage_path:                              ( data.ncbi_lineage_path )
                busco_lineages_folder:                          ( data.busco_lineages_folder )

        }
        .set{ group }


    emit:
    pacbio_reads                     = group.pacbio_reads
    reference                        = group.assembly_path
    assembly_title                   = group.assembly_title
    taxid                            = group.taxid
    nt_database                      = group.nt_database
    nt_kraken_db_path                = group.nt_kraken_db_path
    ncbi_taxonomy_path               = group.ncbi_taxonomy_path
    ncbi_rankedlineage_path          = group.ncbi_rankedlineage_path
    busco_lineages_folder            = group.busco_lineages_folder
    versions                         = ch_versions.ifEmpty(null)
}

def readYAML( yamlfile ) {
    return new Yaml().load( new FileReader( yamlfile.toString() ) )
}