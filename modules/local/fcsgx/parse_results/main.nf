process PARSE_FCSGX_RESULT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9'
        : 'biocontainers/python:3.9'}"

    input:
    tuple val(meta), path(fcs_gx_reports_folder)
    path ncbi_rankedlineage_path

    output:
    tuple val(meta), path("*.csv"), emit: fcsgxresult
    tuple val("${task.process}"), val('parse_fcsgx_result.py'), eval("parse_fcsgx_result.py -v"), topic: versions, emit: versions_parse_fcsgx_result

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    parse_fcsgx_result.py ${fcs_gx_reports_folder} ${ncbi_rankedlineage_path} > ${prefix}_parsed_fcsgx.csv
    """

    stub:
    """
    touch parsed_fcsgx.csv
    """
}
