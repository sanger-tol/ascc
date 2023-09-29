process SED_SED {
    tag "${meta.id}.mummer"
    label "process_medium"

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
            'docker.io/ubuntu:20.04' }"

    input:
    tuple val(meta), path(in_fasta)

    output:
    tuple val(meta), path("*_fixed.fa")   , emit: mummer
    path "versions.yml"                 , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    // args to be '/>/s/ //g'
    """
    sed $args ${in_fasta} > ${prefix}_fixed.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ubuntu: \$(ubuntu --version | sed 's/Ubuntu //g')
        coreutils: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.mummer

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ubuntu: \$(ubuntu --version | sed 's/Ubuntu //g')
        coreutils: $VERSION
    END_VERSIONS
    """
}
