process CREATE_BTK_DATASET {
    label 'process_high'

    input:
    tuple val(meta),        path(reference)
    path gc_content,        stageAs: "?/GC_CONTENT.txt"
    path dot_genome,        stageAs: "?/SORTED.genome"
    path kmers,             stageAs: "?/KMERS_dim_reduction_embeddings_combined.csv"
    path tiara,             stageAs: "?/TIARA.txt"
    path nt_blast,          stageAs: "?/BLAST_with_LINEAGE.csv"
    path mito,              stageAs: "?/MITOCHO.contamination_recommendation"
    path chloro,            stageAs: "?/PLASTID.contamination_recommendation"
    path fcs_adapt,         stageAs: "?/*"
    path fcsgx,             stageAs: "?/FCSGX_parsed.csv"
    path barcode,           stageAs: "?/*"
    path coverage,          stageAs: "?/COVERAGE_AVERAGE.txt"
    path vecscreen,         stageAs: "?/VECSCREEN.vecscreen_contamination"
    path kraken_class,      stageAs: "?/KRAKEN_CLASSIFIED.txt"
    path kraken_report,     stageAs: "?/KRAKEN_REPORT.txt"
    path kraken_lineage,    stageAs: "?/KRAKEN_LINEAGE.txt"

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = args.ext.prefix   ?: "${meta.id}"
    def args    = args.ext.args     ?: ""

    """

    mkdir -p btk_datasets/

        create_btk_dataset.py \\
        ${reference} \\
        btk_datasets/ \\
        ./1/ \\
        ${meta.id} \\
        ${meta.sci_name} \\
        ${meta.taxid} \\
        ${nt_blast} \\
        NA \\
        NA \\
        ${coverage} \\
        $


    """
}