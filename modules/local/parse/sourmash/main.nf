process PARSE_SOURMASH {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c4/c452d87bf8039e112d8e868781331d220b0127e0d8c626650894317139cf49d0/data' :
        'community.wave.seqera.io/library/python_pip_polars:8ef3dfabf82d3525' }"


    input:
    tuple val(meta), path(sourmash_multisearch_results, stageAs: "sourmash_results_*")
    path(assembly_taxa_db_files, stageAs: "taxa_db_*")
    val target_taxa

    output:
    tuple val(meta), path("*.summary.csv")    , emit: multisearch_summary
    tuple val(meta), path("*.non_target.csv") , emit: multisearch_non_target
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sourmash_taxonomy_parser.py \\
        ${args} \\
        -s ${sourmash_multisearch_results.join(' ')} \\
        -a ${assembly_taxa_db_files.join(' ')} \\
        --target_taxa $target_taxa \\
        -o ./
    
    # Rename output files to include the sample prefix
    mv sourmash_results*.summary.csv ${prefix}.multisearch_results.summary.csv
    mv sourmash_results*.non_target.csv ${prefix}.multisearch_results.non_target.csv

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash_taxonomy_parser: \$(sourmash_taxonomy_parser.py --version 2>&1 | cut -d' ' -f2)
    END_VERSIONS
    """
}
