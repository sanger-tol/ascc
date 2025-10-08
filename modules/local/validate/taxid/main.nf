process VALIDATE_TAXID {
    tag "$taxid"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    val(taxid)
    path(ncbi_taxonomy_path)

    output:
    path "versions.yml",        emit: versions

    script:
    """
    find_taxid_in_taxdump.py \\
        $taxid \\
        ${ncbi_taxonomy_path}/nodes.dmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        find_taxid_in_taxdump: \$(find_taxid_in_taxdump.py -v)
    END_VERSIONS
    """

    stub:
    """
    OUTPUT="TAXID FOUND

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        find_taxid_in_taxdump: \$(find_taxid_in_taxdump.py -v)
    END_VERSIONS
    """
}
