process SANGER_TOL_BTK {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta),    path(reference)
    path(samplesheet_csv)
    path blastp,                        stageAs: "blastp.dmnd"
    path blastn
    path blastx
    path tax_dump
    path( "input_pacbio_files/*" )
    val busco_lineages_folder
    val busco_lineages
    val taxon

    output:
    tuple val(meta), path("${meta.id}_btk_out/blobtoolkit/${meta.id}*"),    emit: dataset
    path("${meta.id}_btk_out/blobtoolkit/plots"),                           emit: plots
    path("${meta.id}_btk_out/blobtoolkit/${meta.id}*/summary.json.gz"),     emit: summary_json
    path("${meta.id}_btk_out/busco"),                                       emit: busco_data
    path("${meta.id}_btk_out/multiqc"),                                     emit: multiqc_report
    path("blobtoolkit_pipeline_info"),                                      emit: pipeline_info
    path "versions.yml",                                                    emit: versions

    script:
    def prefix              =   task.ext.prefix         ?:  "${meta.id}"
    def args                =   task.ext.args           ?:  ""
    def profiles            =   task.ext.profiles       ?:  ""
    def get_version         =   task.ext.version_data   ?:  "UNKNOWN - SETTING NOT SET"
    def pipeline_version    =   task.ext.version        ?: "0.6.0"
    // Seems to be an issue where a nested pipeline can't see the files in the same directory
    // Running realpath gets around this but the files copied into the folder are
    // now just wasted space. Should be fixed with using Mahesh's method of nesting but
    // this is proving a bit complicated with BTK

    // outdir should be an arg

    // blastx and blastp use the same database hence the StageAs


    """
    nextflow run sanger-tol/blobtoolkit \\
        -r $pipeline_version \\
        -profile  $profiles \\
        --input "\$(realpath $samplesheet_csv)" \\
        --outdir ${prefix}_btk_out \\
        --fasta "\$(realpath $reference)" \\
        --busco $busco_lineages_folder \\
        --busco_lineages $busco_lineages \\
        --taxon $taxon \\
        --taxdump "\$(realpath $tax_dump)" \\
        --blastp "\$(realpath blastp.dmnd)" \\
        --blastn "\$(realpath $blastn)" \\
        --blastx "\$(realpath $blastx)" \\
        --use_work_dir_as_temp true \\
        --align \\
        $args

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
    def pipeline_version    =   task.ext.version        ?: "draft_assemblies"

    """
    mkdir -p ${meta.id}_btk_out/blobtoolkit/${meta.id}_out
    touch ${meta.id}_btk_out/blobtoolkit/${meta.id}_out/test.json.gz

    mkdir ${meta.id}_btk_out/blobtoolkit/plots
    touch ${meta.id}_btk_out/blobtoolkit/plots/test.png

    mkdir ${meta.id}_btk_out/busco
    touch ${meta.id}_btk_out/busco/test.batch_summary.txt
    touch ${meta.id}_btk_out/busco/test.fasta.txt
    touch ${meta.id}_btk_out/busco/test.json

    mkdir ${meta.id}_btk_out/multiqc
    mkdir ${meta.id}_btk_out/multiqc/multiqc_data
    mkdir ${meta.id}_btk_out/multiqc/multiqc_plots
    touch ${meta.id}_btk_out/multiqc/multiqc_report.html

    mv ${meta.id}_btk_out/pipeline_info blobtoolkit_pipeline_info

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Blobtoolkit: $pipeline_version
        Nextflow: \$(nextflow -v | cut -d " " -f3)
        executor system: $get_version
    END_VERSIONS
    """
}
