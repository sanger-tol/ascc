process ORGANELLE_CONTAMINATION_RECOMMENDATIONS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9'
        : 'biocontainers/python:3.9'}"

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta), path("*contamination_recommendation"), emit: recommendations
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('organelle_contamination_recommendation.py'), eval("organelle_contamination_recommendation.py -v"), topic: versions, emit: versions_organelle_contamination_recommendation

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}-${meta.organelle}"
    """
    organelle_contamination_recommendation.py \
        --input "${files}" \
        ${args} \
        --output ${prefix}.contamination_recommendation
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}-${meta.organelle}"
    """
    touch ${prefix}.bed
    """
}
