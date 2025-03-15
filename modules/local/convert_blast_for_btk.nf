process CONVERT_BLAST_FOR_BTK {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(blast_hits)

    output:
    tuple val(meta), path("*_btk_format.tsv"), emit: btk_format
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    convert_blast_for_btk.py --input ${blast_hits} --output ${prefix}_btk_format.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        convert_blast_for_btk: 1.0.0
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_btk_format.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        convert_blast_for_btk: 1.0.0
    END_VERSIONS
    """
}
