process FILTER_BARCODE {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(barcodes)

    output:
    tuple val(meta), path( "*txt" ) , emit: debarcoded
    path "versions.yml"             , emit: versions

    script:
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    def args    = task.ext.args     ?: ''
    def prefix  = task.ext.prefix   ?: "${meta.id}.debarcoded" // args is for `pacbio_multiplexing_barcodes_check_${meta.barcode}.txt`
    """
    filter_barcode_blast_results.py \\
        --input ${fasta} \\
        --barcodes ${barcodes} \\
        --blast ${blast_data} \\
        --output ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        filter_barcode_blast_results: \$(filter_barcode_blast_results.py -v)
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    def args    = task.ext.args     ?: ''
    def prefix  = task.ext.prefix   ?: "${meta.id}.debarcoded" // args is for `pacbio_multiplexing_barcodes_check_${meta.barcode}.txt`
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        filter_barcode_blast_results: \$(filter_barcode_blast_results.py -v)
    END_VERSIONS
    """
}
