process GET_TARGET_TAXA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), val(taxid)
    path( rankedlineage )
    val( taxonomy_level )

    output:
    tuple val(meta), path( "*.txt" ) , emit: target_taxa
    path "versions.yml"              , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    get_target_taxa_from_taxid.py \\
        --taxid ${taxid} \\
        --rankedlineage ${rankedlineage} \\
        --level ${taxonomy_level} \\
        --output ${prefix}_target_taxa.txt

    # Check if extraction was successful
    if grep -q "TAXID_NOT_FOUND\\|LEVEL_EMPTY" ${prefix}_target_taxa.txt; then
        echo "WARNING: Could not extract target taxa for taxid ${taxid} at level ${taxonomy_level}" >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        get_target_taxa_from_taxid: \$(get_target_taxa_from_taxid.py --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' || echo '1.0.0')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "order:Artiodactyla" > ${prefix}_target_taxa.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        get_target_taxa_from_taxid: 1.0.0
    END_VERSIONS
    """
}
