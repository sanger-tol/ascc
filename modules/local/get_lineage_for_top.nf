process GET_LINEAGE_FOR_TOP {
    tag "${meta.id}"
    label 'process_low'

    conda ""
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    '' :
    '' }"

    input:
    tuple val(meta), path(tophits)
    path( ncbi_taxonomy_path )
    path( ncbi_lineage_path )

    output:
    tuple val(meta), path( "*.tsv" ) , emit: full
    path "versions.yml"              , emit: versions

    script:
    """
    python3 get_lineage_for_top.py ${tophits} ./ ${ncbi_taxonomy_path} ${ncbi_lineage_path} --column_name_prefix nt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        get_lineage_for_top: \$(python3 get_lineage_for_top.py -v)
    END_VERSIONS
    """

    stub:
    """
    touch full_coords.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        get_lineage_for_top: \$(python3 get_lineage_for_top.py -v)
    END_VERSIONS
    """
}
