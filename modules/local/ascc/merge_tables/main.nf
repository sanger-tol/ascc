process ASCC_MERGE_TABLES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta),            path (gc_content,               stageAs: "GC.txt"),
        val(meta_coverage),     path (coverage),
        val(meta_tiara),        path (tiara,                    stageAs: "TIARA.txt"),
        val(meta_kraken),       path (nt_kraken,                stageAs: "LINEAGE.txt"),
        val(meta_ntblast),      path (nt_blast),
        val(meta_kmer),         path (dim_reduction_embeddings),
        val(meta_nrdiamond),    path (nr_diamond),
        val(meta_undiamon),     path (uniprot_diamond,          stageAs: "UP_DIAMOND.tsv"),
        val(meta_btk),          path (btk,                      stageAs: "BTK_summary_table_full.tsv"),
        val(meta_busco_btk),    path (btk_busco),
        val(meta_fcs),          path (fcs_gx,                   stageAs: "FCSGX_parsed.csv")

    output:
    tuple val(meta), path("*_contamination_check_merged_table.csv")         , emit: merged_table
    tuple val(meta), path("*_contamination_check_merged_table_extended.csv"), optional: true, emit: extended_table
    tuple val(meta), path("*_phylum_counts_and_coverage.csv")               , optional: true, emit: phylum_counts
    path "versions.yml",                                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix                      = task.ext.prefix           ?: "${meta.id}"
    def args                        = task.ext.args             ?: ""
    def coverage                    = coverage                  ? "-c ${coverage}"                  : ""
    def tiara                       = tiara                     ? "-t ${tiara}"                     : ""
    def nt_kraken                   = nt_kraken                 ? "-nk ${nt_kraken}"                : ""
    def nt_blast                    = nt_blast                  ? "-nb ${nt_blast}"                 : ""
    def dim_reduction_embeddings    = dim_reduction_embeddings  ? "-dr ${dim_reduction_embeddings}" : ""
    def nr_diamond                  = nr_diamond                ? "-nd ${nr_diamond}"               : ""
    def uniprot_diamond             = uniprot_diamond           ? "-ud ${uniprot_diamond}"          : ""
    def btk                         = btk                       ? "-btk ${btk}"                     : ""
    def btk_busco                   = btk_busco                 ? "-bb ${btk_busco}"                : ""
    def fcs_gx                      = fcs_gx                    ? "-fg ${fcs_gx}"                   : ""

    """
    ascc_merge_tables.py \\
        --gc_cov $gc_content \\
        --sample_name $meta.id \\
        $coverage \\
        $tiara \\
        $nt_kraken \\
        $nt_blast \\
        $dim_reduction_embeddings \\
        $nr_diamond \\
        $uniprot_diamond \\
        $btk \\
        $btk_busco \\
        $fcs_gx \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        ascc_merge_tables: \$(ascc_merge_tables.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix                      = task.ext.prefix           ?:  "${meta.id}"
    """
    touch ${prefix}_contamination_check_merged_table.csv
    touch ${prefix}_contamination_check_merged_table_extended.csv
    touch ${prefix}_phylum_counts_and_coverage.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        ascc_merge_tables: \$(ascc_merge_tables.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

}
