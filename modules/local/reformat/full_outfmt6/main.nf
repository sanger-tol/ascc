process REFORMAT_FULL_OUTFMT6 {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9'
        : 'biocontainers/python:3.9'}"

    input:
    tuple val(meta), path(full_blast)

    output:
    tuple val(meta), path("*outfmt6.tsv"), emit: full
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('reformat_blast_outfmt6.py'), eval("reformat_blast_outfmt6.py -v"), topic: versions, emit: versions_reformat_blast_outfmt6

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    reformat_blast_outfmt6.py ${full_blast} > ${prefix}_outfmt6.tsv
    """

    stub:
    """
    touch full_coords.tsv
    """
}
