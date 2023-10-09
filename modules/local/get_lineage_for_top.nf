process GET_LINEAGE_FOR_TOP {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(tophits)
    path( ncbi_taxonomy_path )
    path( ncbi_lineage_path )

    output:
    tuple val(meta), path( "*.tsv" ) , emit: full
    path "versions.yml"              , emit: versions

    script:
    """
    get_lineage_for_top.py ${tophits} ./ ${ncbi_taxonomy_path} ${ncbi_lineage_path} --column_name_prefix nt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        get_lineage_for_top: \$(get_lineage_for_top.py -v)
    END_VERSIONS
    """

    stub:
    """
    touch full_coords.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        get_lineage_for_top: \$(get_lineage_for_top.py -v)
    END_VERSIONS
    """
}
