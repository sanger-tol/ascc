process SANGER_TOL_BTK {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(reference, stageAs: "REFERENCE.fa")
    tuple val(meta1), path(samplesheet_csv, stageAs: "SAMPLESHEET.csv")
    path blastp, stageAs: "blastp.dmnd"
    path blastn
    path blastx
    path btk_config_file
    path tax_dump
    path btk_yaml, stageAs: "BTK.yaml"
    val taxon
    val gca_accession

    output:
    // path("${meta.id}_btk_out/plots"),                       emit: btk_plots
    path("${meta.id}_btk_out/busco"),                       emit: btk_busco
    path("${meta.id}_btk_out/blobtoolkit"),                 emit: btk_dataset
    path("${meta.id}_btk_out/multiqc"),                     emit: btk_multiqc
    path("blobtoolkit_pipeline_info"),                      emit: btk_pipeline

    script:
    def prefix              =   task.ext.prefix         ?:  "${meta.id}"
    def args                =   task.ext.args           ?:  ""
    def executor            =   task.ext.executor       ?:  ""
    def profiles            =   task.ext.profiles       ?:  ""
    def get_version         =   task.ext.version_data   ?:  "UNKNOWN - SETTING NOT SET"
    def btk_config          =   btk_config_file         ? "-c $btk_config_file"         : ""
    def pipeline_version    =   task.ext.version        ?: "main"
    // YAML used to avoid the use of GCA accession number
    //    https://github.com/sanger-tol/blobtoolkit/issues/77

    // Seems to be an issue where a nested pipeline can't see the files in the same directory
    // Running realpath gets around this but the files copied into the folder are
    // now just wasted space.

    // outdir should be an arg

    """
    $executor 'nextflow run sanger-tol/blobtoolkit \\
        -r $pipeline_version \\
        -profile  $profiles \\
        --input "\$(realpath $samplesheet_csv)" \\
        --outdir ${prefix}_btk_out \\
        --fasta "\$(realpath REFERENCE.fa)" \\
        --yaml "\$(realpath BTK.yaml)" \\
        --accession draft \\
        --taxon $taxon \\
        --taxdump "\$(realpath $tax_dump)" \\
        --blastp "\$(realpath blastp.dmnd)" \\
        --blastn "\$(realpath $blastn)" \\
        --blastx "\$(realpath $blastx)" \\
        $btk_config \\
        $args'

    mv ${prefix}_btk_out/pipeline_info blobtoolkit_pipeline_info

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Blobtoolkit: $pipeline_version
        Nextflow: \$(nextflow -v | cut -d " " -f3)
        executor system: $get_version
    END_VERSIONS
    """

    stub:
    def prefix              =   task.ext.prefix         ?:  "${meta.id}"
    def pipeline_version    =   task.ext.version        ?: "main"

    """
    mkdir -p blobtoolkit/$gca_accession
    touch blobtoolkit/$gca_accession/test.json.gz

    mkdir ${prefix}_btk_out/plots
    touch ${prefix}_btk_out/plots/test.png

    mkdir ${prefix}_btk_out/busco
    touch ${prefix}_btk_out/busco/test.batch_summary.txt
    touch ${prefix}_btk_out/busco/test.fasta.txt
    touch ${prefix}_btk_out/busco/test.json

    mkdir ${prefix}_btk_out/multiqc
    mkdir ${prefix}_btk_out/multiqc/multiqc_data
    mkdir ${prefix}_btk_out/multiqc/multiqc_plots
    touch ${prefix}_btk_out/multiqc/multiqc_report.html

    mv ${prefix}_btk_out/pipeline_info blobtoolkit_pipeline_info

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Blobtoolkit: $pipeline_version
        Nextflow: \$(nextflow -v | cut -d " " -f3)
        executor system: $get_version
    END_VERSIONS
    """
}
