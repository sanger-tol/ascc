process CONVERT_TO_HITS_FILE {
    tag "$meta.id"
    label 'process_low'

    container 'sanger-tol/ascc_btk:3.2.6-c1'

    input:
    tuple val(meta), path(blast_full)

    output:
    tuple val(meta), path("*csv"),      emit: hits_file
    path "versions.yml",                emit: versions

    script:
    def args    =   task.ext.args
    """
    convert_to_hits.py $blast_full $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        convert_to_hits: \$(convert_to_hits.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}