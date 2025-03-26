process CREATE_BTK_DATASET {
    tag "$meta.id"
    label 'process_medium'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CREATE_BTK_DATASET module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    container "docker.io/genomehubs/blobtoolkit:4.3.9"

    input:
    tuple val(id), val(meta), path(reference),
        val(nt_blast_meta), path(nt_blast_file, stageAs: "BLAST_HITS.tsv"),
        val(tiara_meta), path(tiara_file, stageAs: "TIARA.txt"),
        val(kraken2_meta), path(kraken2_file, stageAs: "KRAKEN_REPORT.txt"),
        val(genome_meta), path(dot_genome, stageAs: "SORTED.genome"),
        val(kmers_meta), path(kmers_file, stageAs: "KMERS_dim_reduction_embeddings_combined.csv"),
        val(fcsgx_meta), path(fcsgx_file, stageAs: "FCSGX_parsed.csv"),
        val(nr_full_meta), path(nr_full_file, stageAs: "NUCLEOT_DIAMOND_FULL.tsv"),
        val(un_full_meta), path(un_full_file, stageAs: "UNIPROT_DIAMOND_FULL.tsv"),
        val(mapped_bam_meta), path(mapped_bam_file, stageAs: "MAPPED.bam"),
        val(coverage_meta), path(coverage_file, stageAs: "COVERAGE_AVERAGE.txt"),
        val(kraken1_meta), path(kraken1_file, stageAs: "KRAKEN_CLASSIFIED.txt"),
        val(kraken3_meta), path(kraken3_file, stageAs: "KRAKEN_LINEAGE.txt")

    path ncbi_taxdump,  stageAs: "TAXDUMP"

    output:
    tuple val(meta), path("btk_datasets_CBD"),              emit: btk_datasets
    tuple val(meta), path("btk_summary_table_full.tsv"),    emit: create_summary
    path "versions.yml",                                    emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix          = task.ext.prefix   ?: "${meta.id}"
    def args            = task.ext.args     ?: ""
    def blastn_arg      = nt_blast_file     ? "-bh ${nt_blast_file}"     : ""
    def nt_diamond_arg  = nr_full_file      ? "-nr ${nr_full_file}"      : ""
    def un_diamond_arg  = un_full_file      ? "-ud ${un_full_file}"      : ""
    def kraken_arg      = kraken3_file      ? "-k ${kraken3_file}"       : ""
    def mapped_arg      = mapped_bam_file   ? "-r ${mapped_bam_file}"    : ""
    def tiara_arg       = tiara_file        ? "-t ${tiara_file}"         : ""
    def pca_arg         = kmers_file        ? "-p ${kmers_file}"         : ""
    def fcs_arg         = fcsgx_file        ? "-fc ${fcsgx_file}"        : ""

    """
    mkdir -p btk_datasets_CBD/

    create_btk_dataset.py \\
        -o btk_datasets_CBD \\
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
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        create_btk_dataset: \$(create_btk_dataset.py -v)
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix   ?: "${meta.id}"

    """
    mkdir btk_datasets
    touch btk_datasets/${prefix}.txt
    touch btk_summary_table_full.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        create_btk_dataset: \$(create_btk_dataset.py -v)
    END_VERSIONS
    """
}
