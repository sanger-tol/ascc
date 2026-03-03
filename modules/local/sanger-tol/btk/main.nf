process SANGER_TOL_BTK {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(reference), path(samplesheet_csv)
    path blastp, stageAs: "blastp.dmnd"
    path blastn
    path blastx
    path tax_dump
    path pacbio_files
    path blobtoolkit_config_file
    path blobtoolkit_trace_config
    val busco_lineages_folder
    val busco_lineages
    val taxon

    output:
    tuple val(meta), path("${prefix}_btk_out/blobtoolkit/${prefix}*"),      emit: dataset
    path("${prefix}_btk_out/blobtoolkit/plots"),                            emit: plots
    path("${prefix}_btk_out/blobtoolkit/${prefix}*/summary.json.gz"),       emit: summary_json
    path("${prefix}_btk_out/busco"),                                        emit: busco_data
    path("${prefix}_btk_out/multiqc"),                                      emit: multiqc_report
    path("blobtoolkit_pipeline_info"),                                      emit: pipeline_info
    path "versions.yml",                                                    emit: versions

    script:
    prefix                  =   task.ext.prefix             ?: "${meta.id}"
    def args                =   task.ext.args               ?: ""
    def profiles            =   task.ext.profiles           ?: ""
    def get_version         =   task.ext.version_data       ?: "UNKNOWN - SETTING NOT SET"
    def pipeline            =   task.ext.pipeline           ?: "sanger-tol/blobtoolkit"
    def pipeline_version    =   task.ext.version            ?: "0.8.0"
    def trace_config        =   blobtoolkit_trace_config    ? "-c ${blobtoolkit_trace_config}" : ""
    def singularity_cache_dir = System.getenv('NXF_SINGULARITY_CACHEDIR') ?: System.getenv('SINGULARITY_CACHEDIR') ?: System.getenv('APPTAINER_CACHEDIR') ?: ""
    // Only redirect the image cache here. Do not override TMPDIR:
    // extracting SIF to a sandbox can be very slow on some filesystems, and /tmp is often faster.
    def cache_exports = singularity_cache_dir ? """
    export NXF_SINGULARITY_CACHEDIR="${singularity_cache_dir}"
    export SINGULARITY_CACHEDIR="${singularity_cache_dir}"
    export APPTAINER_CACHEDIR="${singularity_cache_dir}"
    """ : ""
    // Seems to be an issue where a nested pipeline can't see the files in the same directory
    // Running realpath gets around this but the files copied into the folder are
    // now just wasted space. Should be fixed with using Mahesh's method of nesting but
    // this is proving a bit complicated with BTK

    // outdir should be an arg

    // blastx and blastp use the same database hence the StageAs

    // First rename the input fasta, this is done to remove the "_filtered" string
    // from the fasta name.
    // Without doing this, and without the end user knowing, the BTK viewer can
    // end up with non-functional enteries.
    // e.g. uploading the entry iyTipFemo_PRIMARY (which is what the resulting dataset is named)
    // will be a blank entry, with iyTipFemo_PRIMARY_filtered being the correct name due to the input fasta
    // both could be blank because of this.
    //
    """
    mv $reference ${prefix}.fasta

    $cache_exports

    nextflow run ${pipeline} \\
        -r $pipeline_version \\
        -c $blobtoolkit_config_file \\
        ${trace_config} \\
        -profile ${profiles} \\
        --input "\$(realpath $samplesheet_csv)" \\
        --outdir ${prefix}_btk_out \\
        --fasta ${prefix}.fasta \\
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

    # BTK v0.9.0 can name the dataset directory by accession (e.g. GCA_...),
    # while ASCC expects a ${prefix}-based directory for downstream channel matching.
    dataset_dir=\$(find ${prefix}_btk_out/blobtoolkit -mindepth 1 -maxdepth 1 -type d ! -name plots | head -n 1)
    if [[ -n "\$dataset_dir" ]]; then
        ln -sfn "\$(basename "\$dataset_dir")" "${prefix}_btk_out/blobtoolkit/${prefix}"
    fi

    mv ${prefix}_btk_out/pipeline_info blobtoolkit_pipeline_info

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Blobtoolkit: $pipeline_version
        Nextflow: \$(nextflow -v | cut -d " " -f3)
        executor system: $get_version
    END_VERSIONS
    """

    stub:
    prefix              =   task.ext.prefix         ?: "${meta.id}"
    def pipeline_version    =   task.ext.version        ?: "draft_assemblies"

    """
    mkdir -p ${prefix}_btk_out/blobtoolkit/${prefix}_out
    touch ${prefix}_btk_out/blobtoolkit/${prefix}_out/test.json.gz

    mkdir ${prefix}_btk_out/blobtoolkit/plots
    touch ${prefix}_btk_out/blobtoolkit/plots/test.png

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
    END_VERSIONS
    """
}
