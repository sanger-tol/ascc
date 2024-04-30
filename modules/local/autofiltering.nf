process AUTOFILTER_ASSEMBLY {
    tag "$meta.id"
    label "process_medium"

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta),        path(reference)
    tuple val(tiara_meta),  path(tiara_txt)
    tuple val(fcs_meta),    path(fcs_csv)

    output:
    tuple val(meta), path("*autofiltered.fasta"),                       emit: decontaminated_assembly
    tuple val(meta), path("fcs-gx_and_tiara_combined_summary.csv"),     emit: fcs_tiara_summary
    tuple val(meta), path("assembly_filtering_removed_sequences.txt"),  emit: removed_seqs
    path("fcs-gx_alarm_indicator_file.txt"),                            emit: alarm_file

    script:
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    def args    = task.ext.args     ?: ""
    """
    autofilter.py \\
        $reference \\
        $meta.taxid \\
        $tiara_txt \\
        $fcs_csv

    abnormal_contamination_check.py \\
        $reference \\
        fcs-gx_and_tiara_combined_summary.csv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
