include { GENERATE_HTML_REPORT } from '../../modules/local/generate_html_report'
import groovy.json.JsonOutput

workflow GENERATE_HTML_REPORT_WORKFLOW {
    take:
    barcode_results            // channel: [ val(meta), [ barcode_results ] ]
    fcs_adaptor_euk            // channel: [ val(meta), [ fcs_adaptor_euk ] ]
    fcs_adaptor_prok           // channel: [ val(meta), [ fcs_adaptor_prok ] ]
    trim_ns_results            // channel: [ val(meta), [ trim_ns_results ] ]
    vecscreen_results          // channel: [ val(meta), [ vecscreen_results ] ]
    autofilter_results         // channel: [ val(meta), [ autofilter_results ] ]
    merged_table               // channel: [ val(meta), [ merged_table ] ]
    kmers_results              // channel: [ val(meta), [ kmers_results_dirs ] ]
    reference_fasta            // channel: [ val(meta), [ reference_fasta ] ]
    fasta_sanitation_log       // channel: [ val(meta), [ fasta_sanitation_log ] ]
    fasta_length_filtering_log // channel: [ val(meta), [ fasta_length_filtering_log ] ]
    jinja_template             // channel: [ jinja_template ]
    samplesheet                // channel: [ samplesheet ]
    params_file                // channel: [ params_file ]

    main:
    ch_versions = Channel.empty()

    // Add debug logging
    log.info "HTML Report Generation - Input Channels:"
    log.info "barcode_results: ${barcode_results.dump()}"
    log.info "fcs_adaptor_euk: ${fcs_adaptor_euk.dump()}"
    log.info "fcs_adaptor_prok: ${fcs_adaptor_prok.dump()}"
    log.info "trim_ns_results: ${trim_ns_results.dump()}"
    log.info "vecscreen_results: ${vecscreen_results.dump()}"
    log.info "autofilter_results: ${autofilter_results.dump()}"
    log.info "merged_table: ${merged_table.dump()}"
    log.info "fasta_sanitation_log: ${fasta_sanitation_log.dump()}"
    log.info "fasta_length_filtering_log: ${fasta_length_filtering_log.dump()}"
    log.info "jinja_template: ${jinja_template.dump()}"

    // Use mix instead of join to handle empty channels
    // Create a channel with just the meta information
    barcode_results
        .map { meta, files -> meta }
        .unique()
        .set { meta_channel }

    // Map each input channel to a tuple with meta and files, using empty list for missing files
    meta_channel
        .combine(barcode_results.map { meta, files -> [meta.id, files] }.ifEmpty { [null, []] })
        .combine(fcs_adaptor_euk.map { meta, files -> [meta.id, files] }.ifEmpty { [null, []] })
        .combine(fcs_adaptor_prok.map { meta, files -> [meta.id, files] }.ifEmpty { [null, []] })
        .combine(trim_ns_results.map { meta, files -> [meta.id, files] }.ifEmpty { [null, []] })
        .combine(vecscreen_results.map { meta, files -> [meta.id, files] }.ifEmpty { [null, []] })
        .combine(autofilter_results.map { meta, files -> [meta.id, files] }.ifEmpty { [null, []] })
        .combine(merged_table.map { meta, files -> [meta.id, files] }.ifEmpty { [null, []] })
        .combine(kmers_results.map { meta, files -> [meta.id, files] }.ifEmpty { [null, []] })
        .combine(reference_fasta.map { meta, files -> [meta.id, files] }.ifEmpty { [null, []] })
        .combine(fasta_sanitation_log.map { obj -> 
            if (obj instanceof Map) {
                return [obj.id, obj.file ?: []]
            } else if (obj instanceof List && obj.size() >= 2) {
                def meta = obj[0]
                def files = obj[1]
                return [meta.id ?: null, files ?: []]
            } else {
                return [null, []]
            }
        }.ifEmpty { [null, []] })
        .combine(fasta_length_filtering_log.map { obj -> 
            if (obj instanceof Map) {
                return [obj.id, obj.file ?: []]
            } else if (obj instanceof List && obj.size() >= 2) {
                def meta = obj[0]
                def files = obj[1]
                return [meta.id ?: null, files ?: []]
            } else {
                return [null, []]
            }
        }.ifEmpty { [null, []] })
        .map { meta, barcode_id, barcode, fcs_euk_id, fcs_euk, fcs_prok_id, fcs_prok, trim_ns_id, trim_ns, vecscreen_id, vecscreen, autofilter_id, autofilter, merged_id, merged, kmers_id, kmers, ref_id, ref, sanitation_id, sanitation, length_filtering_id, length_filtering ->
            tuple(meta, barcode, fcs_euk, fcs_prok, trim_ns, vecscreen, autofilter, merged, kmers, ref, sanitation, length_filtering)
        }
        .set { combined_inputs }

    log.info "Combined inputs: ${combined_inputs.dump()}"

    // Convert params to JSON for passing to the HTML report
    def paramsJson = JsonOutput.toJson(params)
    
    // Generate HTML report
    GENERATE_HTML_REPORT (
        combined_inputs,
        jinja_template,
        samplesheet,
        params_file,
        paramsJson
    )
    ch_versions = ch_versions.mix(GENERATE_HTML_REPORT.out.versions)

    emit:
    report    = GENERATE_HTML_REPORT.out.report
    versions  = ch_versions
}
