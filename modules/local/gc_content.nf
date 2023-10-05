process GC_CONTENT {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9 conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path( "*-gc.txt" ) , emit: txt
    path "versions.yml"               , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gc_content.py ${fasta} > ${prefix}-gc.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gc_content: \$(gc_content.py -v)
    END_VERSIONS
    """

    stub:
    """
    touch full_coords.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gc_content: \$(gc_content.py -v)
    END_VERSIONS
    """
}
