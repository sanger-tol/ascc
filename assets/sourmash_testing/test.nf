#!/usr/bin/env nextflow

include { RUN_SOURMASH } from '../../subworkflows/local/run_sourmash/main.nf'

workflow {
    // Create input channels
    genome_fasta = Channel.fromPath(params.genome_fasta)
        .map { fasta ->
            def meta = [id: fasta.getBaseName()]
            [meta, fasta]
        }

    sourmash_database = Channel.fromList(params.sourmash_databases)
    assembly_taxa_db = Channel.fromPath(params.assembly_taxa_db)
    target_taxa = Channel.value(params.target_taxa)
    k = params.k
    s = params.s

    // Run the subworkflow
    RUN_SOURMASH(
        genome_fasta,
        sourmash_database,
        assembly_taxa_db,
        target_taxa,
        k,
        s
    )
}
