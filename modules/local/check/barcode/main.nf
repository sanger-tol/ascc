process CHECK_BARCODE {
    tag "${barcodes}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta) , path(pacbio_dir, stageAs:"in/*")
    path barcodes
    val multiplex_csv

    output:
    env OUTPUT          , emit: result
    path "versions.yml" , emit: versions

    script:
    def prefix      = task.ext.prefix   ?: "${meta.id}"
    def args        = task.ext.args     ?: ''
    """
    OUTPUT=\$(\\
        pacbio_barcode_check.py \\
            -b ${barcodes} \\
            -p in/ \\
            -m ${multiplex_csv})

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pacbio_barcode_check: \$(pacbio_barcode_check.py -v)
    END_VERSIONS
    """

    stub:
    """
    OUTPUT="BARCODES FOUND"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pacbio_barcode_check: \$(pacbio_barcode_check.py -v)
    END_VERSIONS
    """
}
