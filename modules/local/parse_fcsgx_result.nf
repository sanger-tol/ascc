process PARSE_FCSGX_RESULT {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(fcs_gx_reports_folder)
    path ncbi_rankedlineage_path

    output:
    tuple val(meta), path( "*.csv" ), emit: fcsgxresult
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    parse_fcsgx_result.py ${fcs_gx_reports_folder} ${ncbi_rankedlineage_path} > parsed_fcsgx.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parse_fcsgx_result: \$(parse_fcsgx_result.py -v)
    END_VERSIONS
    """

    stub:
    """
    touch parsed_fcsgx.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parse_fcsgx_result: \$(parse_fcsgx_result.py -v)
    END_VERSIONS
    """
}
