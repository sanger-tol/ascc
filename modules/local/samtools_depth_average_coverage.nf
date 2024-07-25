process SAMTOOLS_DEPTH_AVERAGE_COVERAGE {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(depth)

    output:
    tuple val(meta), path( "*.txt" ),   emit: average_coverage
    path "versions.yml",                emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools_depth_average_coverage.py $depth > ${prefix}_average_coverage.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        samtools_depth_average_coverage: \$(samtools_depth_average_coverage.py --version)
    END_VERSIONS
    """

    stub:
    """
    touch ${prefix}_average_coverage.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        samtools_depth_average_coverage: \$(samtools_depth_average_coverage.py --version)
    END_VERSIONS
    """
}
