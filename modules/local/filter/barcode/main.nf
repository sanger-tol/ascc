process FILTER_BARCODE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/biopython:1.81'
        : 'biocontainers/biopython:1.81'}"

    input:
    tuple val(meta), path(fasta), path(blast_data), val(barcodes)

    output:
    tuple val(meta), path("*filtered.txt"), emit: debarcoded
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('biopython'), eval("python3 -c 'import Bio; print(Bio.__version__)'"), topic: versions, emit: versions_biopython
    tuple val("${task.process}"), val('filter_barcode_blast_results.py'), eval("filter_barcode_blast_results.py -v"), topic: versions, emit: versions_filter_barcode_blast_results

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    filter_barcode_blast_results.py \
        --input ${fasta} \
        --barcode ${barcodes} \
        --blast ${blast_data} \
        --output ${prefix}_${barcodes}_filtered.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def barcodes_stub = "bc1008"
    """
    touch ${prefix}_${barcodes_stub}_filtered.txt
    """
}
