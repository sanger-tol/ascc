process AUTOFILTER_AND_CHECK_ASSEMBLY {
    tag "${meta.id}"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta),        path(reference)
    tuple val(tiara_meta),  path(tiara_txt)
    tuple val(fcs_meta),    path(fcs_csv)
    path ncbi_rankedlineage

    output:
    tuple val(meta), path("*autofiltered.fasta"),                       emit: decontaminated_assembly
    tuple val(meta), path("*ABNORMAL_CHECK.csv"),                       emit: fcs_tiara_summary
    tuple val(meta), path("assembly_filtering_removed_sequences.txt"),  emit: removed_seqs
    tuple val(meta), path("fcs-gx_alarm_indicator_file.txt"),           emit: alarm_file
    path("*raw_report.txt"),                                            emit: raw_report
    path "versions.yml",                                                emit: versions

    script:
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    def args    = task.ext.args     ?: ""
    def args2   = task.ext.args2    ?: ""
    """
    autofilter.py \\
        ${reference} \\
        --taxid ${meta.taxid} \\
        --tiara ${tiara_txt} \\
        --fcsgx_sum ${fcs_csv} \\
        --out_prefix ${prefix} \\
        --ncbi_rankedlineage_path ${ncbi_rankedlineage} \\
        ${args} \\

    abnormal_contamination_check.py \\
        ${reference} \\
        ${args2} \\
        ${prefix}_ABNORMAL_CHECK.csv \\
        --out_prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch autofiltered.fasta
    touch ABNORMAL_CHECK.csv
    touch assembly_filtering_removed_sequences.txt
    touch fcs-gx_alarm_indicator_file.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
