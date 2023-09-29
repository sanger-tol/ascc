process BLAST_GET_TOP_HITS {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9 conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(outfmt6)

    output:
    tuple val(meta), path( "*tophits.tsv" ) , emit: tophits
    path "versions.yml"                     , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    blast_get_top_hits.py ${outfmt6} > ${prefix}_tophits.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast_get_top_hits: \$(blast_get_top_hits.py -v)
    END_VERSIONS
    """

    stub:
    """
    touch full_coords.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        reformat_blast_outfmt6: \$(blast_get_top_hits.py -v)
    END_VERSIONS
    """
}
