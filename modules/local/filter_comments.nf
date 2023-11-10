process FILTER_COMMENTS {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
            'docker.io/ubuntu:20.04' }"

    input:
    tuple val(meta), path(blast_results)

    output:
    tuple val(meta), path( "*txt" ) , emit: txt
    path "versions.yml"             , emit: versions

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    // EXPLAINER: AWK with space as delimiter, if first value of line is not '#' and column 4 (length of hit) >= 200 then print line and save to file
    """
    awk -v OFS=" " '{ if (\$1 ~ /^[^#]/ && \$4 >= 200) {print \$0} }' ${blast_results} > ${prefix}_filtered_busco.txt

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
    touch ${prefix}.fa

    cat <<-END_VERSIONS > ${prefix}_filtered_busco.txt
    "${task.process}":
        ubuntu: \$(ubuntu --version | sed 's/Ubuntu //g')
        coreutils: $VERSION
    END_VERSIONS
    """
}
