process GC_CONTENT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9'
        : 'biocontainers/python:3.9'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*-GC_CONTENT.txt"), emit: txt
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('gc_content.py'), eval("gc_content.py -v"), topic: versions, emit: versions_gc_content

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gc_content.py ${fasta} > ${prefix}-GC_CONTENT.txt
    """

    stub:
    """
    touch full_coords.tsv
    """
}
