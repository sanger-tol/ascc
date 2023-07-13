process REFORMAT_FULL_OUTFMT6 {
    tag: "${meta.id}"
    label 'process_low'

    conda ""
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    '' :
    '' }"

    input:
    tuple val(meta), path(full_blast)

    output:
    tuple val(meta), path( "*outfmt6.tsv" ) , emit: full
    path "versions.yml"                     , emit: versions

    script:
    """
    python3 reformat_blast_outfmt6.py ${full_blast} > ${meta.id}_outfmt6.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        reformat_blast_outfmt6: \$(python3 reformat_blast_outfmt6.py -v)
    END_VERSIONS
    """

    stub:
    """
    touch full_coords.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        reformat_blast_outfmt6: \$(python3 reformat_blast_outfmt6.py -v)
    END_VERSIONS
    """
}
