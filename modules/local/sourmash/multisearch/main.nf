process SOURMASH_MULTISEARCH {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"

    // Container with sourmash_plugin_branchwater v0.9.14 from Seqera Community
    // See: https://community.wave.seqera.io/
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/81764ecd441b22bc7d6a4694380d710fe71ec8d375a34e2f80de4584cb462c45/data':
        'community.wave.seqera.io/library/sourmash_plugin_branchwater:0.9.14--082047fc8376854f' }"

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
