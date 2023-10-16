process GET_KMER_COUNTS {

    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.9 conda-forge::kcounter=0.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kcounter:0.1.1--py39hf95cd2a_4 ' :
        'quay.io/biocontainers/kcounter:0.0.1' }"

    input:
    tuple val(meta), path(input_fasta)
    val kmer_size

    output:
    path '*_kmer_counts.csv', emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = args.ext.prefix ?: "${meta.id}"
    """
    get_kmers_counts.py \\
        $input_fasta \\
        ${prefix}_kmer_counts.csv \\
        --kmer_size $kmer_size

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        kcounter: \$( kcounter --version 2>&1 | sed 's/kcounter //g' )
        general_purpose_functions.py: \$(general_purpose_functions.py --version | cut -d' ' -f2)
        get_kmers_counts.py: \$(get_kmers_counts.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = args.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_kmer_counts.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        kcounter: \$( kcounter --version 2>&1 | sed 's/kcounter //g' )
        general_purpose_functions.py: \$(general_purpose_functions.py --version | cut -d' ' -f2)
        get_kmers_counts.py: \$(get_kmers_counts.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}



