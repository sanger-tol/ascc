process PARSE_SOURMASH {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d9b9334c4a0777c7722cbcc301a10ddc8684a85f:796f1cea66b7720fa29583d1f6b90404f90dde2f-0' :
        'biocontainers/mulled-v2-d9b9334c4a0777c7722cbcc301a10ddc8684a85f:796f1cea66b7720fa29583d1f6b90404f90dde2f-0' }"


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
