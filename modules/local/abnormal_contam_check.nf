process ABNORMAL_CONTAM_CHECK {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta1), path(fcsgx_tiara_sum)

    output:
    tuple val(meta), path("ABNORMAL_CHECK.csv"),   emit: abnormal

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix              = task.ext.prefix       ?: "${meta.id}"

    """
    abnormal_contamination_check.py \\
        $fasta \\
        $fcsgx_tiara_sum

    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        abnormal_contam_check: \$(abnormal_contam_check.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
