process CHECK_NT_BLAST_TAXONOMY {
    tag "${nt_blast_db_path.baseName}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    path(nt_blast_db_path)

    output:
    path "versions.yml", emit: versions
    stdout emit: status  // Capture the status message from stdout

    script:
    """
    check_nt_blast_taxonomy.py --db_path ${nt_blast_db_path}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        check_nt_blast_taxonomy: \$(check_nt_blast_taxonomy.py --version | sed 's/^//')
    END_VERSIONS
    """

    stub:
    """
    echo "nt_database_taxonomy_files_found"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.0
        check_nt_blast_taxonomy: 1.0.0
    END_VERSIONS
    """
}
