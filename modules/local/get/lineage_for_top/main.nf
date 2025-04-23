process GET_LINEAGE_FOR_TOP {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(tophits)
    path( ncbi_lineage_path )

    output:
    tuple val(meta), path( "*.csv" ) , emit: full
    path "versions.yml"              , emit: versions

    script:
    """
    get_lineage_for_top.py --blast_tsv ${tophits} --taxdump ${ncbi_lineage_path} --output_csv ./${meta.id}_BLAST_results_with_lineage.csv --column_name_prefix nt_

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        get_lineage_for_top: \$(get_lineage_for_top.py --version | sed 's/^//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_BLAST_results_with_lineage.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        get_lineage_for_top: 1.0.0
    END_VERSIONS
    """
}
