process BLAST_CHUNK_TO_FULL {
    tag: "${meta.id}"
    label 'process_low'

    conda ""
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    '' :
    '' }"

    input:
    tuple val(meta), path(chunked)

    output:
    tuple val(meta), path( "*.tsv" ) , emit: full
    path "versions.yml"                             , emit: versions

    script:
    """
    python3 blast_hit_chunk_coords_to_full_coords.py ${chunked}/${meta.id}_nt_blast_chunks.out ${args} > full_coords.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast_hit_chunk_coords_to_full_coords: \$(python3 blast_hit_chunk_coords_to_full_coords.py -v)
    END_VERSIONS
    """

    stub:
    """
    touch full_coords.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast_hit_chunk_coords_to_full_coords: \$(python3 blast_hit_chunk_coords_to_full_coords.py -v)
    END_VERSIONS
    """
}
