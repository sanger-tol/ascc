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
    fcsgx_report_txt           // channel: [ val(meta), [ fcsgx_report_txt ] ]
    fcsgx_taxonomy_rpt         // channel: [ val(meta), [ fcsgx_taxonomy_rpt ] ]

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
    log.info "fcsgx_report_txt: ${fcsgx_report_txt.dump()}"
    log.info "fcsgx_taxonomy_rpt: ${fcsgx_taxonomy_rpt.dump()}"

    // Use mix instead of join to handle empty channels
    // Create a channel keyed by meta.id containing the meta map
    // First, check if barcode_results is empty or has the format [[id: "empty"],[]]
    def meta_channel_by_id
    if (barcode_results.ifEmpty(true)) {
        // If barcode_results is empty, use reference_fasta as the source of meta
        meta_channel_by_id = reference_fasta
            .map { meta, files -> [meta.id, meta] } // Key by meta.id, value is meta map
            .unique()
    } else {
        meta_channel_by_id = barcode_results
            .filter { meta, files -> meta.id != "empty" } // Filter out empty placeholders
            .map { meta, files -> [meta.id, meta] } // Key by meta.id, value is meta map
            .unique()
    }

    // Map each input channel to a tuple keyed by meta.id, using empty list for missing files
    meta_channel_by_id // Now starts with [id, meta]
        .combine(barcode_results.filter { meta, files -> meta.id != "empty" }.map { meta, files -> [meta.id, files] }.ifEmpty { [null, []] }) // Key by meta.id
        .combine(fcs_adaptor_euk.filter { meta, files -> meta.id != "empty" }.map { meta, files -> [meta.id, files] }.ifEmpty { [null, []] }) // Key by meta.id
        .combine(fcs_adaptor_prok.filter { meta, files -> meta.id != "empty" }.map { meta, files -> [meta.id, files] }.ifEmpty { [null, []] }) // Key by meta.id
        .combine(trim_ns_results.filter { meta, files -> meta.id != "empty" }.map { meta, files -> [meta.id, files] }.ifEmpty { [null, []] })
        .combine(vecscreen_results.filter { meta, files -> meta.id != "empty" }.map { meta, files -> [meta.id, files] }.ifEmpty { [null, []] })
        .combine(autofilter_results.filter { meta, files -> meta.id != "empty" }.map { meta, files -> [meta.id, files] }.ifEmpty { [null, []] })
        .combine(merged_table.filter { meta, files -> meta.id != "empty" }.map { meta, files -> [meta.id, files] }.ifEmpty { [null, []] })
        .combine(kmers_results.filter { meta, files -> meta.id != "empty" }.map { meta, files -> [meta.id, files] }.ifEmpty { [null, []] })
        .combine(reference_fasta.map { meta, files -> [meta.id, files] }.ifEmpty { [null, []] }) // Reference fasta should always have a valid meta
        .combine(fasta_sanitation_log.map { obj -> 
            if (obj instanceof Map) {
                return [obj.id, obj.file ?: []]
            } else if (obj instanceof List && obj.size() >= 2) {
                def meta = obj[0]
                // Skip if meta.id is "empty"
                if (meta.id == "empty") {
                    return [null, []]
                }
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
                // Skip if meta.id is "empty"
                if (meta.id == "empty") {
                    return [null, []]
                }
                def files = obj[1]
                return [meta.id ?: null, files ?: []]
            } else {
                return [null, []]
            }
        }.ifEmpty { [null, []] })
        // Handle fcsgx_report_txt and fcsgx_taxonomy_rpt channels which are already in the format [id, file]
        .combine(fcsgx_report_txt.filter { obj -> 
            if (obj instanceof List && obj.size() >= 1) {
                return obj[0] != "empty"
            } else if (obj instanceof Map && obj.containsKey('id')) {
                return obj.id != "empty"
            } else {
                return true
            }
        }.ifEmpty { [null, []] })
        .combine(fcsgx_taxonomy_rpt.filter { obj -> 
            if (obj instanceof List && obj.size() >= 1) {
                return obj[0] != "empty"
            } else if (obj instanceof Map && obj.containsKey('id')) {
                return obj.id != "empty"
            } else {
                return true
            }
        }.ifEmpty { [null, []] })
         // The map closure now receives id, meta, then pairs of id_dup, files for each combined channel
        .map { id, meta, barcode_id_dup, barcode, fcs_euk_id_dup, fcs_euk, fcs_prok_id_dup, fcs_prok, trim_ns_id_dup, trim_ns, vecscreen_id_dup, vecscreen, autofilter_id_dup, autofilter, merged_id_dup, merged, kmers_id_dup, kmers, ref_id_dup, ref, sanitation_id_dup, sanitation, length_filtering_id_dup, length_filtering, fcsgx_report_id_dup, fcsgx_report, fcsgx_tax_id_dup, fcsgx_tax ->
            // Ensure we handle potential nulls from ifEmpty
            def final_barcode = barcode ?: []
            def final_fcs_euk = fcs_euk ?: []
            def final_fcs_prok = fcs_prok ?: []
            def final_trim_ns = trim_ns ?: []
            def final_vecscreen = vecscreen ?: []
            def final_autofilter = autofilter ?: []
            def final_merged = merged ?: []
            def final_kmers = kmers ?: []
            def final_ref = ref ?: []
            def final_sanitation = sanitation ?: []
            def final_length_filtering = length_filtering ?: []
            def final_fcsgx_report = fcsgx_report ?: []
            def final_fcsgx_tax = fcsgx_tax ?: []
            // We only need the original meta map and the actual file data for the output tuple
            tuple(meta, final_barcode, final_fcs_euk, final_fcs_prok, final_trim_ns, final_vecscreen, final_autofilter, final_merged, final_kmers, final_ref, final_sanitation, final_length_filtering, final_fcsgx_report, final_fcsgx_tax)
        }
        .set { combined_inputs }

    log.info "Combined inputs for HTML report: ${combined_inputs.dump()}"

    // Convert params to JSON for passing to the HTML report
    def paramsJson = JsonOutput.toJson(params)
    
    // Generate HTML report
    GENERATE_HTML_REPORT (
        combined_inputs,
        jinja_template,
        samplesheet,
        params_file,
        paramsJson // Note: combined_inputs now contains the fcsgx file paths
    )
    ch_versions = ch_versions.mix(GENERATE_HTML_REPORT.out.versions)

    emit:
    report    = GENERATE_HTML_REPORT.out.report
    versions  = ch_versions
}
