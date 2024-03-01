process KMER_COUNT_DIM_REDUCTION_COMBINE_CSV {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.9 conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(input_files)

    output:
    path '*_kmers_dim_reduction_embeddings_combined.csv', emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = args.ext.prefix ?: "${meta.id}"
    """
    kmer_count_dim_reduction_combine_csv.py \\
        -o ${prefix}_kmers_dim_reduction_embeddings_combined.csv \\
        -i $input_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
        kmer_count_dim_reduction_combine_csv.py: \$(kmer_count_dim_reduction.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = args.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_kmers_dim_reduction_embeddings_combined.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
        kmer_count_dim_reduction.py: \$(kmer_count_dim_reduction.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
