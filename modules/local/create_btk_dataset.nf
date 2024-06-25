process CREATE_BTK_DATASET {
    tag "$meta.id"
    label 'process_medium'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CREATE_BTK_DATASET module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    container "docker.io/genomehubs/blobtoolkit:4.3.9"

    input:
    tuple val(meta),        path(reference)
    path dot_genome,        stageAs: "SORTED.genome"
    path kmers,             stageAs: "KMERS_dim_reduction_embeddings_combined.csv"
    path tiara,             stageAs: "TIARA.txt"
    path nt_blast,          stageAs: "BLAST_HITS.tsv"
    path fcsgx,             stageAs: "FCSGX_parsed.csv"
    path mapped_bam,        stageAs: "MAPPED.bam"
    path coverage,          stageAs: "COVERAGE_AVERAGE.txt"
    path kraken_class,      stageAs: "KRAKEN_CLASSIFIED.txt"
    path kraken_report,     stageAs: "KRAKEN_REPORT.txt"
    path kraken_lineage,    stageAs: "KRAKEN_LINEAGE.txt"
    path nt_diamond,        stageAs: "NUCLEOT_DIAMOND_FULL.tsv"
    path un_diamond,        stageAs: "UNIPROT_DIAMOND_FULL.tsv"
    path ncbi_taxdump,      stageAs: "TAXDUMP"


    /*
            -f ${reference} \\
        -d ./1/ \\
        -n "${prefix}" \\
        -tn "${meta.sci_name}" \\
        -id ${meta.taxid} \\
        -td ${ncbi_taxdump}/ \\
        $blastn_arg \\
        $nt_diamond_arg \\
        $un_diamond_arg \\
        $kraken_arg \\
        $mapped_arg \\
        $tiara_arg \\
        $pca_arg \\
        $fcs_arg \\
        $args */

    output:
    tuple val(meta), path("btk_datasets"),                  emit: btk_datasets
    tuple val(meta), path("btk_summary_table_full.tsv"),    emit: create_summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix          = task.ext.prefix   ?: "${meta.id}"
    def args            = task.ext.args     ?: ""
    def blastn_arg      = nt_blast          ? "-bh ${nt_blast}"     : ""
    def nt_diamond_arg  = nt_diamond        ? "-nr ${nt_diamond}"   : ""
    def un_diamond_arg  = un_diamond        ? "-ud ${un_diamond}"   : ""
    def kraken_arg      = kraken_lineage    ? "-k ${kraken_lineage}": ""
    def mapped_arg      = mapped_bam        ? "-r ${mapped_bam}"    : ""
    def tiara_arg       = tiara             ? "-t ${tiara}"         : ""
    def pca_arg         = kmers             ? "-p ${kmers}"         : ""
    def fcs_arg         = fcsgx             ? "-fc ${fcsgx}"        : ""
    def marker_arg      = ""
    def contigviz_arg   = ""

    """
    mkdir -p btk_datasets/

    create_btk_dataset_V2.py \\
        -f ${reference} \\
        -d ./1/ \\
        -n "${prefix}" \\
        -tn "${meta.sci_name}" \\
        -id ${meta.taxid} \\
        -td ${ncbi_taxdump}/ \\
        $blastn_arg \\
        $nt_diamond_arg \\
        $un_diamond_arg \\
        $kraken_arg \\
        $mapped_arg \\
        $tiara_arg \\
        $pca_arg \\
        $fcs_arg \\
        $args\\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        create_btk_dataset: \$(general_purpose_functions.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
