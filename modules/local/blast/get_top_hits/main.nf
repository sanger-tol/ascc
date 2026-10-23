process BLAST_GET_TOP_HITS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pandas:1.5.2'
        : 'quay.io/biocontainers/pandas:1.5.2'}"

    input:
    tuple val(meta), path(outfmt6)

    output:
    tuple val(meta), path("*tophits.tsv"), emit: tophits
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('v'), eval("blast_get_top_hits.py -v"), topic: versions, emit: versions_blast_get_top_hits

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    blast_get_top_hits.py ${outfmt6} > ${prefix}_tophits.tsv
    """

    stub:
    """
    touch full_coords.tsv
    """
}
