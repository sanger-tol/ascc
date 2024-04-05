process CREATE_BTK_DATASET {
    label 'process_medium'

    input:
    tuple val(meta),        path(reference)
    path dot_genome,        stageAs: "?/SORTED.genome"
    path kmers,             stageAs: "?/KMERS_dim_reduction_embeddings_combined.csv"
    path tiara,             stageAs: "?/TIARA.txt"
    path nt_blast,          stageAs: "?/BLAST_HITS.tsv"
    path fcsgx,             stageAs: "?/FCSGX_parsed.csv"
    path mapped_bam,        stageAs: "?/MAPPED.bam"
    path coverage,          stageAs: "?/COVERAGE_AVERAGE.txt"
    path kraken_class,      stageAs: "?/KRAKEN_CLASSIFIED.txt"
    path kraken_report,     stageAs: "?/KRAKEN_REPORT.txt"
    path kraken_lineage,    stageAs: "?/KRAKEN_LINEAGE.txt"
    path nt_diamond,        stageAs: "?/NUCLEOT_DIAMOND_FULL.tsv"
    path un_diamond,        stageAs: "?/UNIPROT_DIAMOND_FULL.tsv"
    path ncbi_taxdump

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = args.ext.prefix   ?: "${meta.id}"
    def args    = args.ext.args     ?: ""

    """
    mkdir -p btk_datasets/
    ls -lh

    create_btk_dataset.py \\
        ${reference} \\
        btk_datasets/ \\
        ./1/ \\
        ${meta.id} \\
        ${meta.sci_name} \\
        ${meta.taxid} \\
        ${nt_blast} \\
        ${nt_diamond} \\
        ${un_diamond} \\
        ${mapped_bam} \\
        ${ncbi_taxdump}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        create_btk_dataset: \$(general_purpose_functions.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}