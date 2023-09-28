process REFORMAT_FULL_OUTFMT6 {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9 conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(diamond_blast)

    output:
    tuple val(meta), path( "*_diamond_outfmt6.tsv" ) , emit: full
    path "versions.yml"                     , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    reformat_diamond_outfmt6.py ${diamond_blast} > ${prefix}_diamond_outfmt6.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        reformat_diamond_outfmt6: \$(reformat_diamond_outfmt6.py -v)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_diamond_outfmt6.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        reformat_diamond_outfmt6: \$(reformat_diamond_outfmt6.py -v)
    END_VERSIONS
    """
}
