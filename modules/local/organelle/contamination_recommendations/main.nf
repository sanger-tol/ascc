process ORGANELLE_CONTAMINATION_RECOMMENDATIONS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta), path( "*contamination_recommendation" ) , emit: recommendations
    path "versions.yml"                                      , emit: versions

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}-${meta.organelle}"
    """
    organelle_contamination_recommendation.py \\
        --input "${files}" \\
        --output ${prefix}.contamination_recommendation

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        organelle_contamination_recommendation: \$(organelle_contamination_recommendation.py -v)
    END_VERSIONS
    """

    stub:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}-${meta.organelle}"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        organelle_contamination_recommendation: \$(organelle_contamination_recommendation.py -v)
    END_VERSIONS
    """
}
