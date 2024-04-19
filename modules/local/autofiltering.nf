process AUTOFILTER_ASSEMBLY {
    tag "$meta.id"
    label "process_medium"

    container 'docker://quay.io/sanger-tol/ascc_main:0.001-c1'

    input:
    tuple val(meta),        path(reference)
    tuple val(tiara_meta),  path(tiara_txt)
    tuple val(fcs_meta),    path(fcs_csv)

    output:
    tuple val(meta), path("*autofiltered.fasta"),                       emit: decontaminated_assembly
    tuple val(meta), path("fcs-gx_and_tiara_combined_summary.csv"),     emit: fcs_tiara_summary
    tuple val(meta), path("assembly_filtering_removed_sequences.txt")   emit: removed_seqs

    script:
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    def args    = task.ext.args     ?: ""
    """
    remove_fcs_gx_and_tiara.py \\
        $reference \\
        $meta.taxid \\
        $tiara_txt \\
        $fcs_csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}