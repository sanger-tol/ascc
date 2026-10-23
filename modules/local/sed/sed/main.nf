process SED_SED {
    tag "${meta.id}"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
            'docker.io/ubuntu:20.04' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fa")   , emit: sed
    path "versions.yml"             , emit: versions

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def args    = task.ext.args ?: ""
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    sed $args ${fasta} > ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    """
}
