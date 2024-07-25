process GENERATE_SAMPLESHEET {
    tag "$meta.id"
    label "process_low"

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta),    path(pacbio_path)

    output:
    tuple val(meta),    path("*csv"),       emit: csv
    path "versions.yml",                   emit: versions

    script:
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    def args    = task.ext.args     ?: ""
    """
    generate_samplesheet.py \\
        $prefix \\
        $pacbio_path

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        generate_samplesheet: \$(generate_samplesheet.py -v)
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix   ?: "${meta.id}"

    """
    touch ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generate_samplesheet: \$(generate_samplesheet.py -v)
    END_VERSIONS
    """
}
