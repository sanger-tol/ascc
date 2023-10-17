process KMER_COUNT_DIM_REDUCTION {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-82f8313bfec228dfb77c1fa6363fb6be268e81e5:06e3fcd88e927d176149c5bdb78d3a5c7a688be1-0' :
        'biocontainers/mulled-v2-82f8313bfec228dfb77c1fa6363fb6be268e81e5:06e3fcd88e927d176149c5bdb78d3a5c7a688be1-0' }"

    input:
    tuple val(meta), path(kmer_counts_file)
    val dimensionality_reduction_methods
    val n_neighbors_setting
    val autoencoder_epochs_count

    output:
    path '*_kmers_dim_reduction_embeddings.csv', emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = args.ext.prefix ?: "${meta.id}"
    """
    get_lineage_for_kraken_results.py \\
        $kmer_counts_file \\
        ${prefix}_kmers_dim_reduction_embeddings.csv \\
        --selected_methods $ncbi_rankedlineage_path \\
        --n_neighbors_setting $n_neighbors_setting \\
        --autoencoder_epochs_count $autoencoder_epochs_count

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        tensorflow: \$(tensorflow --version | sed 's/tensorflow //g')
        scikit-learn: \$(scikit-learn --version | sed 's/scikit-learn //g')
        umap-learn: \$(umap-learn --version | sed 's/umap-learn //g')
        matplotlib: \$(matplotlib --version | sed 's/matplotlib //g')
        kmer_count_dim_reduction.py: \$(kmer_count_dim_reduction.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = args.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_kmers_dim_reduction_embeddings.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        tensorflow: \$(tensorflow --version | sed 's/tensorflow //g')
        scikit-learn: \$(scikit-learn --version | sed 's/scikit-learn //g')
        umap-learn: \$(umap-learn --version | sed 's/umap-learn //g')
        matplotlib: \$(matplotlib --version | sed 's/matplotlib //g')
        kmer_count_dim_reduction.py: \$(kmer_count_dim_reduction.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
