process SOURMASH_MULTISEARCH {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/81764ecd441b22bc7d6a4694380d710fe71ec8d375a34e2f80de4584cb462c45/data':
        'community.wave.seqera.io/library/sourmash_plugin_branchwater:0.9.14--082047fc8376854f' }"

    input:
    tuple val(meta), path(signature), path(db), val(k), val(s)

    output:
    tuple val(meta), path("*multisearch_results.csv"), emit: multisearch_results
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args   ?: ""
    def prefix         = task.ext.prefix ?: "${meta.id}"
    def PLUGIN_VERSION = "0.9.14" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def out_csv        = "${prefix}_vs_${db.getBaseName()}.multisearch_results.csv"
    """
    sourmash scripts multisearch \\
        ${args} \\
        -k $k \\
        -s $s \\
        -m DNA \\
        -o tmp_multisearch_raw.csv \\
        $signature \\
        $db

    # Keep header + rows where intersect_hashes (column 9) > 0
    awk -F',' 'NR==1 || \$9+0 > 0' tmp_multisearch_raw.csv > "${out_csv}"
    rm tmp_multisearch_raw.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
        sourmash_plugin_branchwater: ${PLUGIN_VERSION}
    END_VERSIONS
    """

    stub:
    def prefix         = task.ext.prefix ?: "${meta.id}"
    def PLUGIN_VERSION = "0.9.14" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.multisearch_results.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
        sourmash_plugin_branchwater: ${PLUGIN_VERSION}
    END_VERSIONS
    """
}
