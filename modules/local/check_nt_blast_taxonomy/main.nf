process CHECK_NT_BLAST_TAXONOMY {
    tag "Check nt BLAST taxonomy"
    label 'process_low'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    path(nt_blast_db_path)

    output:
    path "versions.yml", emit: versions

    script:
    """
    check_nt_blast_taxonomy.py --db_path ${nt_blast_db_path}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        check_nt_blast_taxonomy: \$(check_nt_blast_taxonomy.py --version | sed 's/^//')
    END_VERSIONS
    """
}
