process MERGE_BTK_DATASETS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "docker.io/genomehubs/blobtoolkit:4.3.9"

    input:
    tuple val(meta), path(create_btk_datasets), path(btk_busco_datasets)

    output:
    tuple val(meta), path("merged_btk_datasets"),                                   emit: merged_datasets
    tuple val(meta), path("merged_btk_datasets/btk_busco_summary_table_full.tsv"),  emit: busco_summary_tsv
    path "versions.yml",                                                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix              = task.ext.prefix       ?: "${meta.id}"
    def args                = task.ext.args         ?: ""

    """
    mkdir -p merged_datasets/

    merge_btk_datasets.py \\
        -m $create_btk_datasets \\
        -o ./merged_btk_datasets \\
        -b $btk_busco_datasets \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        merge_btk_datasets: \$(merge_btk_datasets.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix   ?: "${meta.id}"

    """
    mkdir -p merged_datasets/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merge_btk_datasets: \$(merge_btk_datasets.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
