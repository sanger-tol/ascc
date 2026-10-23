process FILTER_COMMENTS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ubuntu:20.04'
        : 'docker.io/ubuntu:20.04'}"

    input:
    tuple val(meta), path(blast_results)

    output:
    tuple val(meta), path("*txt"), emit: txt
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('coreutils'), val("9.1"), topic: versions, emit: versions_coreutils

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    // EXPLAINER: AWK with space as delimiter, if first value of line is not '#' and column 4 (length of hit) >= 200 then print line and save to file
    """
    awk -v OFS=" " '{ if (\$1 ~ /^[^#]/ && \$4 >= 200) {print \$0} }' ${blast_results} > ${prefix}_filtered_blast_result.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa
    """
}
