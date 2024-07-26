process MERGE_BTK_DATASETS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(create_btk_datasets)
    tuple val(meta2), path(busco_btk_datasets)

    output:
    tuple val(meta), path("merged_datasets"),                                   emit: merged_datasets
    tuple val(meta), path("merged_datasets/btk_busco_summary_table_full.tsv"),  emit: busco_summary_tsv
    path "versions.yaml",                                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix              = task.ext.prefix       ?: "${meta.id}"
    def args                = task.ext.args         ?: ""

    """
    mkdir -p merged_datasets/

    merge_btk_datasets_V2.py \\
        -m $create_btk_datasets \\
        -o ./merged_datasets \\
        -b $busco_btk_datasets \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        merge_btk_datasets_V2: \$(merge_btk_datasets_V2.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix   ?: "${meta.id}"

    """
    mkdir -p merged_datasets/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merge_btk_datasets_V2: \$(merge_btk_datasets_V2.py -v)
    END_VERSIONS
    """
}
