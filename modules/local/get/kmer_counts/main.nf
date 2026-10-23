process GET_KMER_COUNTS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/kmer-counter:0.1.2--h4349ce8_0'
        : 'quay.io/biocontainers/kmer-counter:0.1.2--h4349ce8_0'}"

    input:
    tuple val(meta), path(input_fasta)
    val kmer_size

    output:
    tuple val(meta), path(input_fasta), path("*.npy"), emit: npy
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('kmer-counter'), val("0.1.2"), topic: versions, emit: versions_kmer_counter

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    kmer-counter -f ${input_fasta} -k ${kmer_size} -o ${prefix}_kmer_counts.npy
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_kmer_counts.npy
    """
}
