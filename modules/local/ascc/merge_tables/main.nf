process ASCC_MERGE_TABLES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    // rewrite on new input
    tuple val(meta), path(reference),
        path (gc_content),
        path (coverage),
        path (tiara),
        path (nt_kraken),
        path (nt_blast),
        path (dim_reduction_embeddings),
        path (nr_diamond,               stageAs: "diamond.csv"),
        path (uniprot_diamond,          stageAs: "uniprot.csv"),
        path (btk),
        path (btk_busco),
        path (fcs_gx)

    output:
    tuple val(meta), path("*_contamination_check_merged_table.csv")         , emit: merged_table
    tuple val(meta), path("*_contamination_check_merged_table_extended.csv"), optional: true, emit: extended_table
    tuple val(meta), path("*_phylum_counts_and_coverage.csv")               , optional: true, emit: phylum_counts
    path "versions.yml",                                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                             = task.ext.args                    ?: ""

    def empty_file_size                  = 80
    def coverage_data                    = coverage.size()                   > empty_file_size   ? "-c ${coverage}"                  : ""
    def tiara_data                       = tiara.size()                      > empty_file_size   ? "-t ${tiara}"                     : ""
    def nt_kraken_data                   = nt_kraken.size()                  > empty_file_size   ? "-nk ${nt_kraken}"                : ""
    def nt_blast_data                    = nt_blast.size()                   > empty_file_size   ? "-nb ${nt_blast}"                 : ""
    def dim_reduction_embeddings_data    = dim_reduction_embeddings.size()   > empty_file_size   ? "-dr ${dim_reduction_embeddings}" : ""
    def nr_diamond_data                  = nr_diamond.size()                 > empty_file_size   ? "-nd ${nr_diamond}"               : ""
    def uniprot_diamond_data             = uniprot_diamond.size()            > empty_file_size   ? "-ud ${uniprot_diamond}"          : ""
    def btk_data                         = btk.size()                        > empty_file_size   ? "-btk ${btk}"                     : ""
    def btk_busco_data                   = btk_busco.size()                  > empty_file_size   ? "-bb ${btk_busco}"                : ""
    def fcs_gx_data                      = fcs_gx.size()                     > empty_file_size   ? "-fg ${fcs_gx}"                   : ""

    """
    ascc_merge_tables.py \\
        --gc_cov $gc_content \\
        --sample_name $meta.id \\
        $coverage_data \\
        $tiara_data \\
        $nt_kraken_data \\
        $nt_blast_data \\
        $dim_reduction_embeddings_data \\
        $nr_diamond_data \\
        $uniprot_diamond_data \\
        $btk_data \\
        $btk_busco_data \\
        $fcs_gx_data \\
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
