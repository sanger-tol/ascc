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
    tuple val(meta), path ('*_kmers_dim_reduction_embeddings_combined.csv'),    emit: csv
    path "versions.yml",                                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Print debug information
    echo "Meta information: ${meta}"
    echo "Input files: ${input_files}"
    echo "Number of input files: \$(echo ${input_files} | wc -w)"
    echo "Contents of current directory:"
    ls -la
    
    # Check if input files exist
    for file in ${input_files}; do
        echo "Checking file: \$file"
        if [ -f "\$file" ]; then
            echo "  File exists"
            echo "  File size: \$(stat -c%s \$file) bytes"
            echo "  File type: \$(file \$file)"
            echo "  First few lines:"
            head -n 5 "\$file" || echo "  (Failed to read file)"
        else
            echo "  File does not exist"
        fi
    done
    
    # Run the Python script with verbose output
    kmer_count_dim_reduction_combine_csv.py \\
        -o ${prefix}_kmers_dim_reduction_embeddings_combined.csv \\
        -i ${input_files}
    
    # Check the output file
    echo "Output file:"
    if [ -f "${prefix}_kmers_dim_reduction_embeddings_combined.csv" ]; then
        echo "  File exists"
        echo "  File size: \$(stat -c%s ${prefix}_kmers_dim_reduction_embeddings_combined.csv) bytes"
        echo "  First few lines:"
        head -n 5 "${prefix}_kmers_dim_reduction_embeddings_combined.csv" || echo "  (Failed to read file)"
    else
        echo "  File does not exist"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
        kmer_count_dim_reduction_combine_csv.py: \$(python -c 'print("1.0")')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
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
