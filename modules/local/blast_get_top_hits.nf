process BLAST_GET_TOP_HITS {
    tag "${meta.id}"
    label 'process_low'

    conda ""
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    '' :
    '' }"

    input:
    tuple val(meta), path(outfmt6)

    output:
    tuple val(meta), path( "*tophits.tsv" ) , emit: tophits
    path "versions.yml"                     , emit: versions

    script:
    """
    python3 blast_get_top_hits.py ${outfmt6} > ${meta.id}_tophits.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast_get_top_hits: \$(python3 blast_get_top_hits.py -v)
    END_VERSIONS
    """

    stub:
    """
    touch full_coords.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        reformat_blast_outfmt6: \$(python3 blast_get_top_hits.py -v)
    END_VERSIONS
    """
}
