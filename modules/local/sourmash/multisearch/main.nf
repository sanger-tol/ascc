process SOURMASH_MULTISEARCH {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"

    // Custom container with sourmash_plugin_branchwater v0.9.14
    // See: https://hub.docker.com/r/zilov/sourmash_plugin_branchwater
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/zilov/sourmash_plugin_branchwater:0.9.14':
        'docker.io/zilov/sourmash_plugin_branchwater:0.9.14' }"

    input:
    tuple val(meta), path(signature), path(db), val(k), val(s)

    output:
    tuple val(meta), path("*multisearch_results.csv"), emit: multisearch_results
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // required defaults for the tool to run, but can be overridden
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}_vs_${db.getBaseName()}"
    """
    sourmash scripts multisearch \\
        $args \\
        -k $k -s $s -m DNA \\
        -o "${prefix}.multisearch_results.csv" \\
        $signature \\
        $db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_vs_${db.getBaseName()}"
    """
    touch ${prefix}.multisearch_results.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """
}
