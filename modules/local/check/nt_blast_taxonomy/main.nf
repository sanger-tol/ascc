process CHECK_NT_BLAST_TAXONOMY {
    tag "${nt_blast_db_path.baseName}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9'
        : 'biocontainers/python:3.9'}"

    input:
    path nt_blast_db_path

    output:
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('check_nt_blast_taxonomy.py'), eval("check_nt_blast_taxonomy.py --version | sed 's/^//'"), topic: versions, emit: versions_check_nt_blast_taxonomy
    stdout emit: status

    script:
    """
    check_nt_blast_taxonomy.py --db_path ${nt_blast_db_path}
    """

    stub:
    """
    echo "nt_database_taxonomy_files_found"
    """
}
