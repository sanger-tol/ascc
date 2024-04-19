process SANGER_TOL_BTK {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(reference, stageAs: "REFERENCE.fa")
    path samplesheet_csv
    path blastp
    path blastn
    path blastx
    path btk_config
    path tax_dump
    val taxon
    val gca_accession

    output:
    path("blobtoolkit/$gca_accession"), emit: btk_results
    path("blobtoolkit/plots"),          emit: btk_plots
    path("blobktoolkit/busco"),         emit: btk_busco
    path("blobktoolkit/multiqc"),       emit: btk_multiqc
    path("blobtoolkit_pipeline_info"),  emit: btk_pipeline

    script:
    def prefix      =   task.ext.prefix         ?:  "${meta.id}"
    def args        =   task.ext.args           ?:  ""
    def executor    =   task.ext.executor       ?:  ""
    def profiles    =   task.ext.profiles       ?:  ""
    def get_version =   task.ext.version_data   ?:  "UNKNOWN - SETTING NOT SET"
    """
    $executor 'nextflow run sanger-tol/blobtoolkit \\
        -profile  $profiles \\
        --input $samplesheet_csv \\
        --outdir ${meta.id}_btk_out \\
        --fasta $reference \\
        --accession $gc_accession \\
        --taxon $taxon \\
        --taxdump $tax_dump \\
        --blastp $blastp \\
        --blastn $blastn \\
        --blastx $blastx \\
        -c $btk_config \\
        $args'


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Nextflow: \$(nextflow -v | cut -d " " -f3)
        executor system: $get_version
    END_VERSIONS
    """

    stub:
    """
    mkdir -p blobtoolkit/$gca_accession
    touch blobtoolkit/$gca_accession/test.json.gz

    mkdir blobtoolkit/plots
    touch blobtoolkit/plots/test.png

    mkdir blobktoolkit/busco
    touch blobtoolkit/busco/test.batch_summary.txt
    touch blobtoolkit/busco/test.fasta.txt
    touch blobtoolkit/busco/test.json

    mkdir blobktoolkit/multiqc
    mkdir blobktoolkit/multiqc/multiqc_data
    mkdir blobktoolkit/multiqc/multiqc_plots
    touch blobktoolkit/multiqc/multiqc_report.html

    mv pipeline_into blobtoolkit_pipeline_info

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        parse_fcsgx_result: \$(parse_fcsgx_result.py -v)
    END_VERSIONS
    """
}