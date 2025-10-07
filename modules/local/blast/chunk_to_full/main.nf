process BLAST_CHUNK_TO_FULL {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(chunked)

    output:
    tuple val(meta), path( "*.tsv" ) , emit: full
    path "versions.yml"              , emit: versions

    script:
    def args    = task.ext.args     ?: ''
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    """
    $baseDir/bin/blast_hit_chunk_coords_to_full_coords.py ${chunked} ${args} > ${meta.id}_full_coords.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        blast_hit_chunk_coords_to_full_coords: \$(blast_hit_chunk_coords_to_full_coords.py -v)
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix   ?: "${meta.id}"

    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast_hit_chunk_coords_to_full_coords: \$(blast_hit_chunk_coords_to_full_coords.py -v)
    END_VERSIONS
    """
}
