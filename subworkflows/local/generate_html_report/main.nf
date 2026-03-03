include { GENERATE_HTML_REPORT  } from '../../../modules/local/generate_html_report/main'

// FUNCTION IMPORTS
// NOTE: IN FUTURE SHOULD ALSO CONTAIN DATA-MAPPER FUNCTIONS
include { getEmptyPlaceholder   } from '../../../functions/local/ascc_utils'

workflow GENERATE_HTML_REPORT_WORKFLOW {
    take:
    barcode_results            // channel: [ val(meta), [ barcode_results ] ]
    fcs_adaptor
    trim_ns_results            // channel: [ val(meta), [ trim_ns_results ] ]
    vecscreen_results          // channel: [ val(meta), [ vecscreen_results ] ]
    autofilter_results         // channel: [ val(meta), [ autofilter_results ] ]
    merged_table               // channel: [ val(meta), [ merged_table ] ]
    phylum_counts              // channel: [ val(meta), [ phylum_counts ] ] - New input for phylum coverage data
    kmers_results              // channel: [ val(meta), [ kmers_results_dirs ] ]
    reference_fasta            // channel: [ val(meta), [ reference_fasta ] ]
    fasta_sanitation_log       // channel: [ val(meta), [ fasta_sanitation_log ] ]
    fasta_length_filtering_log // channel: [ val(meta), [ fasta_length_filtering_log ] ]
    params_file                // channel: [ params_file ]
    fcsgx_report_txt           // channel: [ val(meta), [ fcsgx_report_txt ] ]
    fcsgx_taxonomy_rpt         // channel: [ val(meta), [ fcsgx_taxonomy_rpt ] ]
    btk_dataset                // channel: [ val(meta), [ btk_dataset ] ]

    main:
    ch_versions = channel.empty()

    // Convert params to JSON for passing to the HTML report
    def paramsJson = groovy.json.JsonOutput.toJson(params)

    //
    // LOGIC: COMBINE ALL INPUT CHANNELS AND ANNOTATE THEM WITH PROCESS TAGS
    //
    reference_fasta
        .map{ meta, file -> [[id: meta.id], file] }

        // REMAINDER means that if there isn't a matching
        // meta.id, it adds a [] in it's place therefore
        // retaining the empty channel
        .join(barcode_results,      remainder: true)
        .join(trim_ns_results,      remainder: true)
        .join(autofilter_results,   remainder: true)
        .join(merged_table,         remainder: true)
        .join(phylum_counts,        remainder: true)
        .join(fcs_adaptor
            .map { meta, files ->
                        // LOGIC: SORT INTO PREDICTABLE ORDER
                        def files_snapshot = files instanceof List ? files.toList() : [files]
                        def sorted_files = files_snapshot.sort { file ->
                            file.toString().contains('_euk') ? 0 :
                            file.toString().contains('_prok') ? 1 : 2
                        }
                        [[id: meta.id], sorted_files]
                    },              remainder: true)
        .join(vecscreen_results,    remainder: true)
        .join(kmers_results,        remainder: true)
        .join(fasta_sanitation_log, remainder: true)
        .join(fasta_length_filtering_log, remainder: true)
        .join(fcsgx_report_txt,     remainder: true)
        .join(fcsgx_taxonomy_rpt,   remainder: true)
        .join(btk_dataset,          remainder: true)

        // FILTER CHAIN IS NEEDED SO THAT WE CAN GET RID OF THE MULTIPLE
        // TYPES OF ERROR CHANNELS WHICH WILL GET OUTPUT FROM THIS
        // EMPTY INPUT CHANNELS CREATE NEW EMPTY OUTPUT CHANNELS
        .filter { items ->
            def meta = items[0]
            meta != null &&
            meta != [] &&
            !(meta instanceof Map && (meta.id == null || meta.isEmpty()))
        }
        // ON THE LEFT OVER CHANNELS, THOSE WITH VALID META OBJECTS
        .map { items ->
            // Replace null values with placeholder file
            items.withIndex().collect { item, index ->
                if (item == null) {
                    getEmptyPlaceholder(index)
                } else if (item instanceof List && item.isEmpty()) {
                    getEmptyPlaceholder(index)
                } else {
                    item
                }
            }
        }
        .set { sorted_data }


    //
    // MODULE: GENERATE A HTML REPORT FOR END USER
    //
    GENERATE_HTML_REPORT (
        sorted_data,
        channel.fromPath("${projectDir}/assets/templates/*.jinja").collect(),  // Pass the list of Jinja templates
        channel.fromPath(params.input).collect(),
        [],   // Do we even use a params.params_file?
        channel.value(paramsJson),                                             // JSON string can be used multiple times
        channel.fromPath("${projectDir}/assets/css/*.css").collect()           // channel.of one (CSS files)
    )
    ch_versions = ch_versions.mix(GENERATE_HTML_REPORT.out.versions)

    emit:
    report    = GENERATE_HTML_REPORT.out.report
    versions  = ch_versions
}
