process REFORMAT_TO_HITS_FILE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9'
        : 'biocontainers/python:3.9'}"

    input:
    tuple val(meta), path(blast_full)

    output:
    tuple val(meta), path("*csv"), emit: hits_file
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('convert_to_hits.py'), eval("convert_to_hits.py --version | cut -d' ' -f2"), topic: versions, emit: versions_convert_to_hits

    script:
    def args = task.ext.args ?: ""
    """
    convert_to_hits.py ${blast_full} ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.csv
    """
}
