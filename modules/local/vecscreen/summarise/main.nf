process SUMMARISE_VECSCREEN_OUTPUT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(vecscreen_filtered_outfile)

    output:
    tuple val(meta), path("*.vecscreen_contamination"), emit: vecscreen_contamination
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    summarise_vecscreen_output.py ${vecscreen_filtered_outfile} > ${prefix}.vecscreen_contamination

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        summarise_vecscreen_output.py: \$(summarise_vecscreen_output.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vecscreen_contamination

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        summarise_vecscreen_output.py: \$(summarise_vecscreen_output.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
