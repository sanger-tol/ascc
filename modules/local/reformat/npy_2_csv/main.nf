process REFORMAT_NPY_2_CSV {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d9b9334c4a0777c7722cbcc301a10ddc8684a85f:796f1cea66b7720fa29583d1f6b90404f90dde2f-0' :
        'biocontainers/mulled-v2-d9b9334c4a0777c7722cbcc301a10ddc8684a85f:796f1cea66b7720fa29583d1f6b90404f90dde2f-0' }"

    input:
    tuple val(meta), path(fasta), path(npy)
    val kmer_size

    output:
    tuple val(meta), path( "*_KMER_COUNTS.csv" ) , emit: csv
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    npy_2_csv.py -f $fasta -k $kmer_size -n $npy -o ${prefix}_KMER_COUNTS.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        npy_2_csv.py: \$(npy_2_csv.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_KMER_COUNTS.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        npy_2_csv.py: \$(npy_2_csv.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}