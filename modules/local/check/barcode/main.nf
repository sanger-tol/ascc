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
    path "barcode_results.txt"  , emit: result
    path "versions.yml"         , emit: versions

    script:
    """
    pacbio_barcode_check.py \\
        -b ${barcodes} \\
        -p in/ \\
        -m ${multiplex_csv} > barcode_results.txt

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
