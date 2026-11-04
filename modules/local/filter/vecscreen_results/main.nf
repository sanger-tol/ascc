process FILTER_VECSCREEN_RESULTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9'
        : 'biocontainers/python:3.9'}"

    input:
    tuple val(meta), path(vecscreen_outfile)

    output:
    tuple val(meta), path("*vecscreen.grepped.out"), emit: filtered_vecscreen_outfile
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('VSlistTo1HitPerLine.py'), eval("VSlistTo1HitPerLine.py -v"), topic: versions, emit: versions_VSlistTo1HitPerLine

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    VSlistTo1HitPerLine.py ${args} ${vecscreen_outfile} > ${prefix}_vecscreen.grepped.out
    """

    stub:
    def prefix = args.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_vecscreen.grepped.out
    """
}
