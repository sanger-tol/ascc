process KMER_COUNT_DIM_REDUCTION_COMBINE_CSV {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pandas:1.5.2'
        : 'quay.io/biocontainers/pandas:1.5.2'}"

    input:
    tuple val(meta), path(input_files)

    output:
    tuple val(meta), path('*_kmers_dim_reduction_embeddings_combined.csv'), emit: csv
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pandas'), eval("python3 -c 'import pandas; print(pandas.__version__)'"), topic: versions, emit: versions_pandas
    tuple val("${task.process}"), val('kmer_count_dim_reduction_combine_csv.py'), eval("kmer_count_dim_reduction.py --version | cut -d' ' -f2"), topic: versions, emit: versions_kmer_count_dim_reduction_combine_csv

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = args.ext.prefix ?: "${meta.id}"
    """
    kmer_count_dim_reduction_combine_csv.py \
        -o ${prefix}_kmers_dim_reduction_embeddings_combined.csv \
        -i ${input_files}
    """

    stub:
    def prefix = args.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_kmers_dim_reduction_embeddings_combined.csv
    """
}
