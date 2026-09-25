process PARSE_SOURMASH {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c4/c452d87bf8039e112d8e868781331d220b0127e0d8c626650894317139cf49d0/data' :
        'community.wave.seqera.io/library/python_pip_polars:8ef3dfabf82d3525' }"


    input:
    // original file names are kept (sub-directory) so the parser can label hits with their DB
    // contig_sizes (name<TAB>length, e.g. the .genome file) or []: every contig gets a classification row
    tuple val(meta), path(sourmash_multisearch_results, stageAs: "sourmash_results/*"), val(target_taxa), path(contig_sizes)
    path(assembly_taxa_db_files, stageAs: "taxa_db_*")
    path(exclude_accessions_file)             // optional: TSV/CSV with assembly_accession column, or []
    path(taxdump, stageAs: "taxdump")         // optional: NCBI new_taxdump dir, for sourmash_div / --exclude_species
    path(div_mapping)                         // optional: BLAST-div -> FCS-GX div table (assets/fcs_gx_blast_div_mapping.tsv)

    output:
    tuple val(meta), path("*.summary.csv")       , emit: multisearch_summary
    tuple val(meta), path("*.non_target.csv")    , emit: multisearch_non_target
    tuple val(meta), path("*.classification.csv"), emit: multisearch_classification
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def exclude_file = exclude_accessions_file ? "--exclude_accessions_file ${exclude_accessions_file}" : ""
    def sizes_arg    = contig_sizes ? "--contig_sizes ${contig_sizes}" : ""
    def taxdump_arg  = taxdump ? "--taxdump ${taxdump}" : ""
    def div_arg      = div_mapping ? "--div_mapping ${div_mapping}" : ""
    """
    sourmash_taxonomy_parser.py \\
        ${args} \\
        ${exclude_file} \\
        ${sizes_arg} \\
        ${taxdump_arg} \\
        ${div_arg} \\
        -s ${sourmash_multisearch_results} \\
        -a ${assembly_taxa_db_files.join(' ')} \\
        --target_taxa $target_taxa \\
        --prefix sourmash_results \\
        -o ./

    # Rename output files to include the sample prefix
    mv sourmash_results*.summary.csv ${prefix}.multisearch_results.summary.csv
    mv sourmash_results*.non_target.csv ${prefix}.multisearch_results.non_target.csv
    mv sourmash_results*.classification.csv ${prefix}.multisearch_results.classification.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash_taxonomy_parser: \$(sourmash_taxonomy_parser.py --version 2>&1 | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.multisearch_results.summary.csv
    touch ${prefix}.multisearch_results.non_target.csv
    touch ${prefix}.multisearch_results.classification.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash_taxonomy_parser: \$(sourmash_taxonomy_parser.py --version 2>&1 | cut -d' ' -f2)
    END_VERSIONS
    """
}
