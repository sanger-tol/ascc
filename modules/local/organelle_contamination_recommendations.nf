process ORGANELLE_CONTAMINATION_RECOMMENDATIONS {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path( "*bed" ) , emit: bed
    path "versions.yml"             , emit: versions

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: ''
    """
    organelle_contamination_recommendation.py \\
        --input ${bed} \\
        --output ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        organelle_contamination_recommendation: \$(organelle_contamination_recommendation.py -v)
    END_VERSIONS
    """

    stub:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: ''
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        organelle_contamination_recommendation: \$(organelle_contamination_recommendation.py -v)
    END_VERSIONS
    """
}
