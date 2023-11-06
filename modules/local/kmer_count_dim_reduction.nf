process KMER_COUNT_DIM_REDUCTION {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ac95cc1cb32439236d915b38af3e056ce8eb0375:2dc02878e5657bb1f787431c323c9d261fc6d520-0' :
        'biocontainers/mulled-v2-ac95cc1cb32439236d915b38af3e056ce8eb0375:2dc02878e5657bb1f787431c323c9d261fc6d520-0' }"

    input:
    tuple val(meta), path(kmer_counts_file)
    val dimensionality_reduction_method
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

    kmer_count_dim_reduction.py \\
        $kmer_counts_file \\
        ${prefix}_${dimensionality_reduction_method}_kmers_dim_reduction_embeddings.csv \\
        --selected_methods $dimensionality_reduction_method \\
        --n_neighbors_setting $n_neighbors_setting \\
        --autoencoder_epochs_count $autoencoder_epochs_count

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
        tensorflow: \$(python3 -c 'import tensorflow; print(tensorflow.__version__)')
        scikit-learn: \$(python3 -c 'import sklearn; sklearn.show_versions()')
        umap-learn: \$(python3 -c 'import umap; print(umap.__version__)')
        matplotlib: \$(python3 -c 'import matplotlib; print(matplotlib.__version__)')
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
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
        tensorflow: \$(tensorflow --version | sed 's/tensorflow //g')
        scikit-learn: \$(scikit-learn --version | sed 's/scikit-learn //g')
        umap-learn: \$(umap-learn --version | sed 's/umap-learn //g')
        matplotlib: \$(matplotlib --version | sed 's/matplotlib //g')
        kmer_count_dim_reduction.py: \$(kmer_count_dim_reduction.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}