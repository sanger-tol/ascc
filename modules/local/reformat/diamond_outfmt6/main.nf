process REFORMAT_DIAMOND_OUTFMT6 {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9'
        : 'biocontainers/python:3.9'}"

    input:
    tuple val(meta), path(diamond_blast)

    output:
    tuple val(meta), path("*_diamond_outfmt6.tsv"), emit: full
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('reformat_diamond_outfmt6.py'), eval("reformat_diamond_outfmt6.py -v"), topic: versions, emit: versions_reformat_diamond_outfmt6

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    reformat_diamond_outfmt6.py ${diamond_blast} > ${prefix}_diamond_outfmt6.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_diamond_outfmt6.tsv
    """
}
