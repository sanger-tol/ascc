process BLAST_CHUNK_TO_FULL {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9'
        : 'biocontainers/python:3.9'}"

    input:
    tuple val(meta), path(chunked)

    output:
    tuple val(meta), path("*_full_coords.tsv"), emit: full
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('blast_hit_chunk_coords_to_full_coords.py'), eval("blast_hit_chunk_coords_to_full_coords.py -v"), topic: versions, emit: versions_blast_hit_chunk_coords_to_full_coords

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    blast_hit_chunk_coords_to_full_coords.py \\
        ${chunked} \\
        ${args} \\
        > ${prefix}_full_coords.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_full_coords.tsv
    """
}
