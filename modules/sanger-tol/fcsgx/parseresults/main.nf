process FCSGX_PARSERESULTS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.14' :
        'quay.io/biocontainers/python:3.14' }"

    input:
    tuple val(meta), path(taxonomy_report)
    tuple val(meta2), path(fcs_report)
    path ncbi_rankedlineage_path

    output:
    tuple val(meta), path( "*_parsed_fcsgx.csv" ),   emit: fcsgxresult
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('parse_fcsgx_result'), eval("parse_fcsgx_result.py --version"), emit: versions_parsefcsgx, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    parse_fcsgx_result.py \\
        ${taxonomy_report} \\
        ${fcs_report} \\
        ${ncbi_rankedlineage_path} > ${prefix}_parsed_fcsgx.csv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_parsed_fcsgx.csv
    """
}
