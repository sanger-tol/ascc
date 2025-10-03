process REFORMAT_TO_HITS_FILE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(blast_full)

    output:
    tuple val(meta), path("*csv"),      emit: hits_file
    path "versions.yml",                emit: versions

    script:
    def args    =   task.ext.args   ?: ""
    """
    $baseDir/bin/convert_to_hits.py $blast_full $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        convert_to_hits: \$(convert_to_hits.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix   ?: "${meta.id}"

    """
    touch ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        convert_to_hits: \$(convert_to_hits.py -v)
    END_VERSIONS
    """
}
