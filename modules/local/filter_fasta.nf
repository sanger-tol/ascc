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
    tuple val(meta), path("fasta_sanitation.json"), emit: sanitation_log
    tuple val(meta), path("fasta_length_filtering.json"), emit: length_filtering_log
    path "versions.yml",                        emit: versions


    script:
    def prefix      = task.ext.prefix   ?: "${meta.id}"
    def args        = task.ext.args     ?: ''
    def max_length  = task.ext.cutoff   ?: 1900000000 // This is the default value and maximum supported length of sequence per scaffold
    """
    sanitise_input_fasta_file.py \\
        ${input_fasta} \\
        --log_file fasta_sanitation.json \\
        --max_detailed_changes 0 > ${prefix}_shortened.fasta

    filter_fasta_by_length.py \\
        ${args} \\
        ${prefix}_shortened.fasta \\
        ${max_length} \\
        --log_file fasta_length_filtering.json > ${prefix}_filtered.fasta 2> ${prefix}_filtering_stderr.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        sanitise_input_fasta_file: \$(sanitise_input_fasta_file.py -v)
        filter_fasta_by_length: \$(filter_fasta_by_length.py -v)
    END_VERSIONS
    """

    stub:
    def prefix      = task.ext.prefix   ?: "${meta.id}"
    """
    touch ${prefix}_filtered.fasta
    touch ${prefix}_filtering_stderr.txt
    echo '{"input_file":"${input_fasta}","total_sequences":0,"has_issues":false}' > fasta_sanitation.json
    echo '{"input_file":"${input_fasta}","total_sequences":0,"sequences_retained":0,"sequences_filtered":0,"has_filtering":false}' > fasta_length_filtering.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        sanitise_input_fasta_file: \$(sanitise_input_fasta_file.py -v)
        filter_fasta_by_length: \$(filter_fasta_by_length.py -v)
    END_VERSIONS
    """
}
