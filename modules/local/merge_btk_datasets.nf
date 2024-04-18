rocess MERGE_BTK_DATASETS {
    tag "$meta.id"
    label 'process_low'

    container 'sanger-tol/ascc_btk:3.2.6-c1'


    input:
    tuple val(meta), path(create_btk_datasets)
    tuple val(meta), path(busco_btk_datasets)
    tuple val(meta), path(busco_summary_file)

    output:
    tuple val(meta), path("merged_datasets"),   emit: merged_datasets

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix              = args.ext.prefix   ?: "${meta.id}"
    def args                = args.ext.args     ?: ""
    def busco_btk_datasets  = busco_btk_datasets? "-b ${busco_btk_datasets}"     : ""
    def nt_diamond_arg      = nt_diamond        ? "-nr ${nt_diamond}"   : ""
    def un_diamond_arg      = un_diamond        ? "-ud ${un_diamond}"   : ""

    """
    mkdir -p merged_datasets/
    ls -lh

    merge_btk_datasets_V2.py \\
        -m $create_btk_datasets \\
        -o ./merged_datasets \\
        $busco_btk_datasets \\
        // $busco_summary_file \\
        $args

    echo "merge_btk_dataset.py btk_datasets/ BTK_FOLDER_PATH BTK_WITH_BUSCO"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        create_btk_dataset: \$(general_purpose_functions.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}