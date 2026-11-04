process VALIDATE_TAXID {
    tag "${taxid}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9'
        : 'biocontainers/python:3.9'}"

    input:
    val taxid
    path ncbi_taxonomy_path

    output:
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('find_taxid_in_taxdump.py'), eval("find_taxid_in_taxdump.py -v"), topic: versions, emit: versions_find_taxid_in_taxdump

    script:
    """
    find_taxid_in_taxdump.py \
        ${taxid} \
        ${ncbi_taxonomy_path}/nodes.dmp
    """

    stub:
    """
    OUTPUT="TAXID FOUND
    """
}
