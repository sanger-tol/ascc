process GET_KMER_COUNTS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'docker.io/ubuntu:20.04' }"

    input:
    tuple val(meta), path(fasta)
    val kmer_size

    output:
    tuple val(meta), path(fasta), path( "*.npy" ), emit: npy
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def KMERCOUNTER_VERSION = "0.1.2"
    def prefix              = task.ext.prefix   ?: "${meta.id}"
    def args                = task.ext.args     ?: ""
    """
    kmer-counter -f $fasta -k $kmer_size -o ${prefix}_output.npy

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmercounter: $KMERCOUNTER_VERSION
    END_VERSIONS
    """

    stub:
    def KCOUNTER_VERSION = "0.1.1"
    def prefix = args.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_KMER_COUNTS.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmercounter: $KMERCOUNTER_VERSION
    END_VERSIONS
    """
}
