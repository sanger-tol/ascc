include { GENERATE_HTML_REPORT } from '../../../modules/local/generate_html_report/main'
import groovy.json.JsonOutput

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
    jinja_templates_list       // channel: [ jinja_templates_list ] // Updated to accept a list
    samplesheet                // channel: [ samplesheet ]
    params_file                // channel: [ params_file ]
    fcsgx_report_txt           // channel: [ val(meta), [ fcsgx_report_txt ] ]
    fcsgx_taxonomy_rpt         // channel: [ val(meta), [ fcsgx_taxonomy_rpt ] ]
    btk_dataset                // channel: [ val(meta), [ btk_dataset ] ]
    css_files_list             // channel: [ css_files_list ] // CSS files to include in the report

    main:
    ch_versions = Channel.empty()

    fcs_adaptor
        .multiMap { meta, file1, file2 ->
            euk: [meta + [process: "FCS-Adaptor-EUK"], file1]
            prok: [meta + [process: "FCS-Adaptor-PROK"], file2]
        }
        .set { fcs_adaptor_split }

    barcode_results
        .map{ meta, _file ->
            def new_meta = meta + [process: "BARCODES"]
            [new_meta, _file]
        }.mix(
            fcs_adaptor_split.euk,
            fcs_adaptor_split.prok,
            trim_ns_results
                .map{ meta, _file ->
                    def new_meta = meta + [process: "TRAILINGNS"]
                    [new_meta, _file]
                },
            vecscreen_results,
            autofilter_results
                .map{ meta, _file ->
                    def new_meta = meta + [process: "AUTOFILTER"]
                    [new_meta, _file]
                },
            merged_table
                .map{ meta, _file ->
                    def new_meta = meta + [process: "MERGED_TABLE"]
                    [new_meta, _file]
                },
            phylum_counts
                .map{ meta, _file ->
                    def new_meta = meta + [process: "MERGED_PHYLUM_COUNTS"]
                    [new_meta, _file]
                },
            kmers_results,
            reference_fasta,
            fasta_sanitation_log,
            fasta_length_filtering_log,
            fcsgx_report_txt,
            fcsgx_taxonomy_rpt,
            btk_dataset
        )
        .map { meta, file ->
            [meta.id, [meta: meta, file: file]]
        }
        .filter { id, data -> id != [] }
        .groupTuple()
        .map { id, data ->
            [id: id, data: data]
        }
        .set { all_data }

        def processes = [
            'REFERENCE', 'REFERENCE_FILT_LOG', 'REFERENCE_SANI_LOG',
            'TRAILING_NS', 'FCSGX_REPORT', 'FCSGX_TAX_REPORT',
            'Vecscreen', 'KMER_RESULTS', "FCS-Adaptor-EUK", "FCS-Adaptor-PROK"
        ]

        def processChannels = processes.collectEntries { process ->
            [(process): all_data
                .map { sample ->
                    def data = sample.data.find { it.meta.process == process }
                    data ? [sample.id, data.meta, data.file] : [sample.id, [process: process], []]
                }
            ]
        }

        def combined_channels = processChannels['REFERENCE']

        processes.tail().each { process ->
            combined_channels = combined_channels
                                    .combine(processChannels[process], by: 0)
        }

        // samplesheet,
        // params_file,
        // jinja_templates_list,
        // css_files_list

        combined_channels.view{"DATA: $it"}

    // Convert params to JSON for passing to the HTML report
    //def paramsJson = JsonOutput.toJson(params)

    // // Generate HTML report
    // GENERATE_HTML_REPORT (
    //     combined_inputs,
    //     jinja_templates_list, // Pass the list of Jinja templates
    //     samplesheet,
    //     params_file,
    //     paramsJson, // Note: combined_inputs now contains the fcsgx file paths
    //     css_files_list // Pass the CSS files
    // )
    // ch_versions = ch_versions.mix(GENERATE_HTML_REPORT.out.versions)

    // emit:
    // report    = GENERATE_HTML_REPORT.out.report
    // versions  = ch_versions
}
