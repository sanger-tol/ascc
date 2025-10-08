process GENERATE_HTML_REPORT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ab48c38c3be93a696d7773767d9287b4a0d3bf19:e3c8a1ac0a27058d7922e8b6d02f303c30d93e3a-0' :
        'quay.io/biocontainers/mulled-v2-ab48c38c3be93a696d7773767d9287b4a0d3bf19:e3c8a1ac0a27058d7922e8b6d02f303c30d93e3a-0' }"

    input:
    tuple val(meta),
        path(barcode_results,               stageAs: "barcodes/*"),
        path(fcs_adaptor_euk,               stageAs: "fcs/*"),
        path(fcs_adaptor_prok,              stageAs: "fcs/*"),
        path(trim_ns_results,               stageAs: "trailingns/*"),
        path(vecscreen_results,             stageAs: "vecscreen/*"),
        path(autofilter_results,            stageAs: "autofilter/*"),
        path(merged_table,                  stageAs: "merged/*"),
        path(phylum_counts,                 stageAs: "coverage/*"),
        path(kmers_results,                 stageAs: "kmers/**"),
        path(reference_fasta),
        path(fasta_sanitation_log,          stageAs: "fasta_sanitation/*"),
        path(fasta_length_filtering_log,    stageAs: "fasta_length_filtering/*"),
        path(fcs_gx_report_txt,             stageAs: "fcsgx/*"),
        path(fcs_gx_taxonomy_rpt,           stageAs: "fcsgx/*"),
        path(btk_output_dir,                stageAs: "btk/*")
    path(jinja_templates_list,              stageAs: "templates/*")
    path(samplesheet)
    path(params_file)
    val(params_json)
    path(css_files_list, stageAs: "css/*")

    output:
    tuple val(meta), path("report/site/*"), emit: report
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix          = task.ext.prefix   ?: "${meta.id}"
    // Convert params_json to a properly escaped string for command line
    def params_json_arg = params_json       ? "--params_json '${params_json.replaceAll("'", "\\'")}'" : ""
    def params_file     = params_file       ? "--params_file ${params_file}" : ""
    def btk_included    = "${params.run_create_btk_dataset == 'both' || (params.run_create_btk_dataset == 'genomic' && params.genomic_only) || (params.run_create_btk_dataset == 'organellar' && !params.genomic_only)}"
    def btk_outdir      = "${params.outdir}/${meta.id}"
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
    generate_html_report.py \\
        --output_dir report \\
        --template_dir report/site/templates \\
        --barcode_dir barcodes \\
        --fcs_dir fcs \\
        --trim_ns_dir trailingns \\
        --vecscreen_dir vecscreen \\
        --autofilter_dir autofilter \\
        --merged_dir merged \\
        --coverage_dir coverage \\
        --kmers_dir kmers \\
        --reference ${reference_fasta} \\
        --fasta_sanitation_log fasta_sanitation/fasta_sanitation.json \\
        --fasta_length_filtering_log fasta_length_filtering/fasta_length_filtering.json \\
        --samplesheet ${samplesheet} \\
        ${params_file} \\
        ${params_json_arg} \\
        --fcs_gx_report_txt fcsgx/*.fcs_gx_report.txt \\
        --fcs_gx_taxonomy_rpt fcsgx/*.taxonomy.rpt \\
        --btk_published_path "${btk_outdir}/create_btk_dataset/btk_datasets_CBD" \\
        --btk_included "${btk_included}" \\
        --launch_dir "${workflow.launchDir}" \\
        --pipeline_version ${workflow.manifest.version} \\
        --output_prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        jinja2: \$(python -c 'import jinja2; print(jinja2.__version__)')
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
        generate_html_report: \$(generate_html_report.py --version 2>&1 | sed 's/.*version //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch report/site/${prefix}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        jinja2: \$(python -c 'import jinja2; print(jinja2.__version__)')
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
        generate_html_report: \$(python generate_html_report.py --version 2>&1 | sed 's/.*version //g')
    END_VERSIONS
    """
}
