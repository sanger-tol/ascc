process MERGE_BTK_DATASETS {
    tag "${meta.id}"
    label 'process_low'
    container "docker.io/genomehubs/blobtoolkit:4.3.9"

    input:
    tuple val(meta), path(create_btk_datasets), path(btk_busco_datasets)

    output:
    tuple val(meta), path("merged_btk_datasets"), emit: merged_datasets
    tuple val(meta), path("merged_btk_datasets/btk_busco_summary_table_full.tsv"), emit: busco_summary_tsv
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('merge_btk_datasets.py'), eval("merge_btk_datasets.py -v"), topic: versions, emit: versions_merge_btk_datasets

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit(1, "MERGE_BTK_DATASETS module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    def args = task.ext.args ?: ""
    """
    mkdir -p merged_datasets/

    merge_btk_datasets.py \\
        -m ${create_btk_datasets} \\
        -o ./merged_btk_datasets \\
        -b ${btk_busco_datasets} \\
        ${args}
    """

    stub:
    """
    mkdir -p merged_datasets/
    """
}
