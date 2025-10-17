process GENERATE_HTML_REPORT {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.outdir}/${meta.id}/html_report", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ab48c38c3be93a696d7773767d9287b4a0d3bf19:e3c8a1ac0a27058d7922e8b6d02f303c30d93e3a-0' :
        'quay.io/biocontainers/mulled-v2-ab48c38c3be93a696d7773767d9287b4a0d3bf19:e3c8a1ac0a27058d7922e8b6d02f303c30d93e3a-0' }"

    input:
    tuple val(meta),
        path(barcode_results, stageAs: "barcodes/*"),
        path(fcs_adaptor_euk, stageAs: "fcs/*"),
        path(fcs_adaptor_prok, stageAs: "fcs/*"),
        path(trim_ns_results, stageAs: "trailingns/*"),
        path(vecscreen_results, stageAs: "vecscreen/*"),
        path(autofilter_results, stageAs: "autofilter/*"),
        path(merged_table, stageAs: "merged/*"),
        path(phylum_counts, stageAs: "coverage/*"),       // Add phylum coverage data input
        path(kmers_results, stageAs: "kmers/**"),
        path(reference_fasta),
        path(fasta_sanitation_log, stageAs: "fasta_sanitation/*"),
        path(fasta_length_filtering_log, stageAs: "fasta_length_filtering/*"),
        path(fcs_gx_report_txt, stageAs: "fcsgx/*"),      // Add FCS-GX report txt input
        path(fcs_gx_taxonomy_rpt, stageAs: "fcsgx/*"),    // Add FCS-GX taxonomy rpt input
        path(btk_output_dir, stageAs: "btk/*")            // Add BTK output dir input
    path(jinja_templates_list, stageAs: "templates/*") // Updated to accept a list and stage them
    path(samplesheet)
    path(params_file)
    val(params_json)
    path(css_files_list, stageAs: "css/*") // CSS files to include in the report

    output:
    tuple val(meta), path("report/site/*.html"), emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Build CLI args using staged file paths (readers). Include only if present.
    def reference_file         = reference_fasta             ? "--reference ${reference_fasta}" : ""
    def barcode_arg            = barcode_results             ? "--barcode_file ${barcode_results}" : ""
    def fcs_euk_arg            = fcs_adaptor_euk             ? "--fcs_adaptor_euk_file ${fcs_adaptor_euk}" : ""
    def fcs_prok_arg           = fcs_adaptor_prok            ? "--fcs_adaptor_prok_file ${fcs_adaptor_prok}" : ""
    def trim_ns_arg            = trim_ns_results             ? "--trim_ns_file ${trim_ns_results}" : ""
    def vecscreen_arg          = vecscreen_results           ? "--vecscreen_file ${vecscreen_results}" : ""
    def autofilter_arg         = autofilter_results          ? "--autofilter_file ${autofilter_results}" : ""
    def merged_table_arg       = merged_table                ? "--merged_table_file ${merged_table}" : ""
    def phylum_coverage_arg    = phylum_counts               ? "--phylum_coverage_file ${phylum_counts}" : ""
    def fasta_sanitation_arg   = fasta_sanitation_log        ? "--fasta_sanitation_log ${fasta_sanitation_log}" : ""
    def fasta_length_filtering_arg = fasta_length_filtering_log ? "--fasta_length_filtering_log ${fasta_length_filtering_log}" : ""
    def samplesheet_arg        = samplesheet                 ? "--samplesheet ${samplesheet}" : ""
    def params_file_arg        = params_file                 ? "--params_file ${params_file}" : ""
    def fcs_gx_report_arg      = fcs_gx_report_txt           ? "--fcs_gx_report_txt ${fcs_gx_report_txt}" : ""
    def fcs_gx_taxonomy_arg    = fcs_gx_taxonomy_rpt         ? "--fcs_gx_taxonomy_rpt ${fcs_gx_taxonomy_rpt}" : ""

    // Convert params_json to a properly escaped string for command line
    def params_json_arg = params_json ? "--params_json '${params_json.replaceAll("'", "\\'")}'" : ""
    """
    mkdir -p report/site/templates
    mkdir -p report/site/css
    cp templates/*.jinja report/site/templates/
    cp css/*.css report/site/css/

    # Create kmers directory if it doesn't exist
    mkdir -p kmers

    # Debug: List the contents of the kmers directory
    echo "Contents of kmers directory:"
    ls -la kmers/

    # Run the report generation script directly from bin directory
    python $baseDir/bin/generate_html_report.py \\
        --output_dir report \\
        --template_dir report/site/templates \\
        $barcode_arg \\
        $fcs_euk_arg \\
        $fcs_prok_arg \\
        $trim_ns_arg \\
        $vecscreen_arg \\
        $autofilter_arg \\
        $merged_table_arg \\
        $phylum_coverage_arg \\
        --kmers_dir kmers \\
        $reference_file \\
        $fasta_sanitation_arg \\
        $fasta_length_filtering_arg \\
        $samplesheet_arg \\
        $params_file_arg \\
        $params_json_arg \\
        $fcs_gx_report_arg \\
        $fcs_gx_taxonomy_arg \\
        --btk_published_path "${params.outdir}/${meta.id}/create_btk_dataset/btk_datasets_CBD" \\
        --btk_included "${params.run_create_btk_dataset == 'both' || (params.run_create_btk_dataset == 'genomic' && params.genomic_only) || (params.run_create_btk_dataset == 'organellar' && !params.genomic_only)}" \\
        --launch_dir "${workflow.launchDir}" \\
        --outdir "${params.outdir}" \\
        --pipeline_version ${workflow.manifest.version} \\
        --output_prefix $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        jinja2: \$(python -c 'import jinja2; print(jinja2.__version__)')
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
        generate_html_report: \$(python $baseDir/bin/generate_html_report.py --version 2>&1 | sed 's/.*version //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p report/site/templates
    mkdir -p report/site/css
    mkdir -p report/site/kmers
    cp -r kmers/* report/site/kmers
    cp css/*.css report/site/css/
    touch report/site/${prefix}.html


    # Run the report generation script (stub mode)
    echo "Would run python $baseDir/bin/generate_html_report.py --pipeline_version ${workflow.manifest.version} --fcs_gx_report_txt fcsgx/*.fcs_gx_report.txt --fcs_gx_taxonomy_rpt fcsgx/*.taxonomy.rpt with samplesheet $samplesheet ${params_file ? 'and params file ' + params_file : ''} and FASTA sanitation logs"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        jinja2: \$(python -c 'import jinja2; print(jinja2.__version__)')
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
        generate_html_report: \$(python $baseDir/bin/generate_html_report.py --version 2>&1 | sed 's/.*version //g')
    END_VERSIONS
    """
}
