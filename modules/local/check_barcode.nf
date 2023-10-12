process CHECK_BARCODE {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta)     , path(barcodes)
    tuple val(meta2)    , path(pacbio_dir)
    tuple val(meta3)    , path(multiplex_csv)

    output:
    stdout              , emit: debarcoded
    path "versions.yml" , emit: versions

    script:
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    def args    = task.ext.args     ?: ''
    """
    pacbio_barcode_check.py \\
        ${barcode_fasta} \\
        ${pacbio_dir} \\
        ${multiplex_csv}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pacbio_barcode_check: \$(pacbio_barcode_check.py -v)
    END_VERSIONS
    """

    stub:
    """
    echo "BARCODES FOUND!"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pacbio_barcode_check: \$(pacbio_barcode_check.py -v)
    END_VERSIONS
    """
}
