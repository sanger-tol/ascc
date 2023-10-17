process FILTER_BARCODE {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9 conda-forge::biopython=1.78"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'biocontainers/biopython:1.78' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(blast_data)
    val barcodes

    output:
    tuple val(meta), path( "*filtered.txt" )    , emit: debarcoded
    path "versions.yml"                         , emit: versions

    script:
    def args    = task.ext.args     ?: ''
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    """
    filter_barcode_blast_results.py \\
        --input ${fasta} \\
        --barcode ${barcodes} \\
        --blast ${blast_data} \\
        --output ${prefix}_${barcodes}_filtered.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python3 -c 'import Bio; print(Bio.__version__)')
        filter_barcode_blast_results: \$(filter_barcode_blast_results.py -v)
    END_VERSIONS
    """

    stub:
    def args        = task.ext.args     ?: ''
    def prefix      = task.ext.prefix   ?: "${meta.id}"
    def barcodes    = "bc1008"
    """
    touch ${prefix}_${barcodes}_filtered.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        filter_barcode_blast_results: \$(filter_barcode_blast_results.py -v)
    END_VERSIONS
    """
}
