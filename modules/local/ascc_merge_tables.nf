process ASCC_MERGE_TABLES {
    tag "$meta.id"
    label 'process_low'

    container 'sanger-tol/ascc_btk:3.2.6-c1'

    input:
    tuple val(meta), path(gc_content, stageAs: "GC.txt")
    path coverage
    path tiara,     stageAs: "TIARA.txt"
    path bacterial_kraken
    path nt_kraken, stageAs: "LINEAGE.txt"
    path nt_blast
    path dim_reduction_embeddings
    path nr_diamond
    path uniprot_diamond,   stageAs: "UP_DIAMOND.tsv"
    path cobiontid_markerscan
    path contigviz
    path btk,        stageAs: "BTK_summary_table_full.tsv"
    path btk_busco
    path fcs_gx,    stageAs: "FCSGX_parsed.csv"

    output:
    tuple val(meta), path("*_contamination_check_merged_table.csv")         , emit: merged_table
    tuple val(meta), path("*_contamination_check_merged_table_extended.csv"), optional: true, emit: extended_table
    tuple val(meta), path("*_phylum_counts_and_coverage.csv")               , optional: true, emit: phylum_counts

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix                      = task.ext.prefix           ?: "${meta.id}"
    def args                        = task.ext.args             ?: ""
    def coverage                    = coverage                  ? "-c ${coverage}"                  : ""
    def tiara                       = tiara                     ? "-t ${tiara}"                     : ""
    def bacterial_kraken            = bacterial_kraken          ? "-bk ${bacterial_kraken}"         : ""
    def nt_kraken                   = nt_kraken                 ? "-nk ${nt_kraken}"                : ""
    def nt_blast                    = nt_blast                  ? "-nb ${nt_blast}"                 : ""
    def dim_reduction_embeddings    = dim_reduction_embeddings  ? "-dr ${dim_reduction_embeddings}" : ""
    def nr_diamond                  = nr_diamond                ? "-nd ${nr_diamond}"               : ""
    def uniprot_diamond             = uniprot_diamond           ? "-ud ${uniprot_diamond}"          : ""
    def contigviz                   = contigviz                 ? "-cv ${contigviz}"                : ""
    def btk                         = btk                       ? "-btk ${btk}"                     : ""
    def btk_busco                   = btk_busco                 ? "-bb ${btk_busco}"                : ""
    def fcs_gx                      = fcs_gx                    ? "-fg ${fcs_gx}"                   : ""
    def cobiontid_markerscan        = ""

    """

    ascc_m_tables.py \\
        --gc_cov $gc_content \\
        --sample_name $meta.id \\
        $coverage \\
        $tiara \\
        $bacterial_kraken \\
        $nt_kraken \\
        $nt_blast \\
        $dim_reduction_embeddings \\
        $nr_diamond \\
        $uniprot_diamond \\
        $contigviz \\
        $btk \\
        $btk_busco \\
        $fcs_gx \\
        $cobiontid_markerscan \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        ascc_merge_tables: \$(ascc_merge_tables.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
