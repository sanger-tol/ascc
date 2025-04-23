process GET_KMER_COUNTS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.9 conda-forge::kcounter=0.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d9b9334c4a0777c7722cbcc301a10ddc8684a85f:796f1cea66b7720fa29583d1f6b90404f90dde2f-0' :
        'biocontainers/mulled-v2-d9b9334c4a0777c7722cbcc301a10ddc8684a85f:796f1cea66b7720fa29583d1f6b90404f90dde2f-0' }"

    input:
    tuple val(meta), path(input_fasta)
    val kmer_size

    output:
    tuple val(meta), path( "*_KMER_COUNTS.csv" ) , emit: csv
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def KCOUNTER_VERSION = "0.1.1"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    get_kmers_counts.py \\
        $input_fasta \\
        ${prefix}_KMER_COUNTS.csv \\
        --kmer_size $kmer_size

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        kcounter: $KCOUNTER_VERSION
        general_purpose_functions.py: \$(general_purpose_functions.py --version | cut -d' ' -f2)
        get_kmers_counts.py: \$(get_kmers_counts.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def KCOUNTER_VERSION = "0.1.1"
    def prefix = args.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_KMER_COUNTS.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        kcounter: $KCOUNTER_VERSION
        general_purpose_functions.py: \$(general_purpose_functions.py --version | cut -d' ' -f2)
        get_kmers_counts.py: \$(get_kmers_counts.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
