process DECONTAMINATE_GENERATE_BED {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81':
        'quay.io/biocontainers/biopython:1.81' }"

    input:
    tuple val(meta),    path(fasta),
        path(parsed_fcsgx,              stageAs: "./dataset/fcsgx_data/*"),
        path(abnormal_check_csv,        stageAs: "./dataset/autofilter/*"),
        path(fcs_adaptor_file,          stageAs: "./dataset/fcs_adaptor/*"),
        path(trailingns,                stageAs: "./dataset/trailingns/*"),
        path(barcodes,                  stageAs: "./dataset/filter_barcode/**"),
        path(sourmash_non_target,       stageAs: "./dataset/sourmash/*"),
        path(mito_recommendations,      stageAs: "./dataset/organelle_contamination_recommendations.mito_mito/*"),
        path(plastid_recommendations,   stageAs: "./dataset/organelle_contamination_recommendations.plastid_plastid/*")
    tuple val(meta2), path(longreads)

    output:
    tuple val(meta), path(fasta), path("*.contamination.bed")   , emit: contamination_bed
    tuple val(meta), path("*.tiara.bed")                        , emit: tiara_bed
    tuple val(meta), path("*.report.txt")                       , emit: jira_report          // Prettified report for JIRA
    tuple val(meta), path("*.abnormal_details.txt")             , emit: abnormal_report
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix      = task.ext.prefix                       ?: "${meta.id}"
    def args        = task.ext.args                         ?: ""
    """
    generate_contamination_bed.py \\
        --assembly_path $fasta \\
        --decon_tolid_type_dir "./dataset" \\
        --assembly_type ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        generate_contamination_bed.py: \$(generate_contamination_bed.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix      = task.ext.prefix   ?: "${meta.id}"
    """
    touch ${prefix}.contamination.bed
    touch ${prefix}.tiara.bed
    touch ${prefix}.report

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        generate_contamination_bed.py: \$(generate_contamination_bed.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
