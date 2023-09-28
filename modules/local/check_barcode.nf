process CHECK_BARCODE {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9 conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"

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
    //barcode_fasta = file path
    //pacbio_dir = /lustre/scratch123/tol/resources/treeval/treeval-testdata/asccTinyTest/pacbio/
    //multiplex_name = "bca101,bcc202"
    """
    pacbio_barcode_check.py \\
        ${barcode_fasta} \\
        ${pacbio_dir} \\
        ${multiplex_csv}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pacbio_barcode_check: \$(pacbio_barcode_check.py -v)
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    def args    = task.ext.args     ?: ''
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pacbio_barcode_check: \$(pacbio_barcode_check.py -v)
    END_VERSIONS
    """
}
