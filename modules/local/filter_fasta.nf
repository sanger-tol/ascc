process FILTER_FASTA {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(input_fasta)

    output:
    tuple val(meta), path("*_filtered.fasta"),  emit: fasta
    path "*_filtering_stderr.txt",              emit: error_log
    path "versions.yml",                                emit: versions


    script:
    def prefix      = task.ext.prefix   ?: "${meta.id}"
    def args        = task.ext.args     ?: ''
    def max_length  = task.ext.cutoff   ?: 1900000000 // This is the default value and maximum supported length of sequence per scaffold
    """
    ascc_shorten_fasta_headers.py \\
        ${input_fasta} > ${prefix}_shortened.fasta

    filter_fasta_by_length.py \\
        ${args} \\
        ${prefix}_shortened.fasta \\
        ${max_length} > ${prefix}_filtered.fasta 2> ${prefix}_filtering_stderr.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        ascc_shorten_fasta_headers: \$(ascc_shorten_fasta_headers.py -v)
        filter_fasta_by_length: \$(filter_fasta_by_length.py -v)
    END_VERSIONS
    """

    stub:
    def prefix      = task.ext.prefix   ?: "${meta.id}"
    """
    touch ${prefix}_filtered.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        ascc_shorten_fasta_headers: \$(ascc_shorten_fasta_headers.py -v)
        filter_fasta_by_length: \$(filter_fasta_by_length.py -v)
    END_VERSIONS
    """
}
