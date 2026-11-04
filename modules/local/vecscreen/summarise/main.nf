process SUMMARISE_VECSCREEN_OUTPUT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9'
        : 'biocontainers/python:3.9'}"

    input:
    tuple val(meta), path(vecscreen_filtered_outfile)

    output:
    tuple val(meta), path("*.vecscreen_contamination"), emit: vecscreen_contamination
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('summarise_vecscreen_output.py'), eval("summarise_vecscreen_output.py --version | cut -d' ' -f2"), topic: versions, emit: versions_summarise_vecscreen_output

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    summarise_vecscreen_output.py ${vecscreen_filtered_outfile} > ${prefix}.vecscreen_contamination
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vecscreen_contamination
    """
}
