process KMER_COUNT_DIM_REDUCTION {

    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-ac95cc1cb32439236d915b38af3e056ce8eb0375:34bd58763b84a3ea2f6c60b87b4858b2f80c070e-0'
        : 'biocontainers/mulled-v2-ac95cc1cb32439236d915b38af3e056ce8eb0375:34bd58763b84a3ea2f6c60b87b4858b2f80c070e-0'}"

    input:
    tuple val(meta), path(kmer_counts_file)
    val dimensionality_reduction_method
    val n_neighbors_setting
    val autoencoder_epochs_count

    output:
    tuple val(meta), path('*_kmers_dim_reduction_dir'), emit: kmers_dim_reduction_dir
    tuple val(meta), path('./'), emit: results_dir
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pandas'), eval("python3 -c 'import pandas; print(pandas.__version__)'"), topic: versions, emit: versions_pandas
    tuple val("${task.process}"), val('tensorflow'), eval("python3 -c 'import tensorflow; print(tensorflow.__version__)'"), topic: versions, emit: versions_tensorflow
    tuple val("${task.process}"), val('scikit-learn'), eval("python3 -c \"import sklearn; print(sklearn.__version__)\""), topic: versions, emit: versions_scikit_learn
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('umap-learn'), val("0.5.5"), topic: versions, emit: versions_umap_learn
    tuple val("${task.process}"), val('matplotlib'), eval("python3 -c 'import matplotlib; print(matplotlib.__version__)'"), topic: versions, emit: versions_matplotlib
    tuple val("${task.process}"), val('kmer_count_dim_reduction.py'), eval("kmer_count_dim_reduction.py --version | cut -d' ' -f2"), topic: versions, emit: versions_kmer_count_dim_reduction

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = args.ext.prefix ?: "${meta.id}"
    def dir_name = "${prefix}_${dimensionality_reduction_method}_kmers_dim_reduction_dir"
    """
    # Create the directory with the standardized naming convention
    mkdir -p ${dir_name}

    # Run the dimensionality reduction and save the output in the directory
    kmer_count_dim_reduction.py \
        ${kmer_counts_file} \
        ${dir_name} \
        --selected_methods ${dimensionality_reduction_method} \
        --n_neighbors_setting ${n_neighbors_setting} \
        --autoencoder_epochs_count ${autoencoder_epochs_count}
    """

    stub:
    def prefix = args.ext.prefix ?: "${meta.id}"
    def dir_name = "${prefix}_${dimensionality_reduction_method}_kmers_dim_reduction_dir"
    """
    mkdir -p ${dir_name}
    touch ${dir_name}/kmers_dim_reduction_embeddings.csv
    """
}
