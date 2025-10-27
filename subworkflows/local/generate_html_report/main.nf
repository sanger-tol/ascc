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

    //
    // LOGIC: COMBINE ALL INPUT CHANNELS AND ANNOTATE THEM WITH PROCESS TAGS
    //
    fcs_adaptor
        .multiMap { meta, file1, file2 ->
            euk: [meta + [process: "FCS_ADAPTOR_EUK"], file1]
            prok: [meta + [process: "FCS_ADAPTOR_PROK"], file2]
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
        // MAP BY THE SAMPLE.ID SO WE HAVE A CLEAR 'KEY'
        .map { meta, file ->
            [meta.id, [meta: meta, file: file]]
        }
        // FILTER AWAY THE NULL CHANNEL
        // THIS GETS MADE NOW THAT ALL CHANNELS HAVE PROCESS TAGS
        // IT FORCES THE CREATION OF THE NULL CHANNEL
        .filter { id, data -> id != [] && id != null }
        // GROUP ON 'KEY' AND RE-MAP INTO A MORE SENSIBLE STRUCTURE
        .groupTuple()
        .map { id, data ->
            [id: id, data: data]
        }
        .set { all_data }

        // LIST OF EXPECTED INPUT DATA
        def processes = [
            'REFERENCE',
            'BARCODES', 'REFERENCE_FILT_LOG', 'REFERENCE_SANI_LOG',
            'TRAILING_NS', 'FCSGX_REPORT', 'FCSGX_TAX_REPORT',
            'VECSCREEN', 'KMER_RESULTS', "FCS_ADAPTOR_EUK",
            "FCS_ADAPTOR_PROK", "AUTOFILTER", "MERGED_TABLE",
            "MERGED_PHYLUM_COUNTS", "BTK_DATASET"
        ]

        // COLLECT THE DATA AND MAP SO THAT CHANNELS NOT EXISTING ARE RECREATED
        // THIS WAY WE WILL ALWAYS HAVE A KNOWN LENGTH CHANNEL, IT'LL BE
        // processes * 2 (TO ACCOUNT FOR META AND FILE)
        def processChannels = processes.collectEntries { process ->
            [(process): all_data
                .map { sample ->
                    def data = sample.data.find { it.meta.process == process }
                    data ? [sample.id, data.meta, data.file] : [sample.id, [process: process], []]
                }
            ]
        }

        // SET THE KEY CHANNEL AND COMBINE EVERYTHING ELSE ON IT
        def combined_channels = processChannels['REFERENCE']
        processes.tail().each { process ->
            combined_channels = combined_channels
                                    .combine(processChannels[process], by: 0)
        }

        // ABOVE PROCESSES ARE NOT DETEMINISTIC, HOWEVER, WE SHOULD BE ABLE TO
        // SORT THEM INTO THE GIVEN ORDER
        combined_channels
            .map { sample ->
                def id = sample[0]  // First element is the common key
                def dataItems = sample[1..-1]  // Rest are the actual data items

                // CREATE A MAP ON THE PROCESSES LIST, THIS WILL BE TRUTH ORDER
                def processOrder = processes.withIndex().collectEntries { proc, idx ->
                    [proc, idx]
                }

                // SORT DATA [meta, file] BY THE PROCESS ORDER
                def sortedData = dataItems.sort { item, file ->

                    def processName = null
                    if (item instanceof List && item.size() > 0 && item[0] instanceof Map) {
                        processName = item[0].process
                    } else if (item instanceof Map) {
                        processName = item.process
                    }
                    // UNKNOWN PROCESSES
                    processOrder[processName] ?: 999
                }

                [[ id: id ], sortedData]
            }
            .map { item ->
                // FORCE INTO [[meta0], meta1, file1, meta2, file2...
                [item[0]] + item[1]
            }
            .set { sorted_data }


    // Convert params to JSON for passing to the HTML report
    def paramsJson = JsonOutput.toJson(params)

    // Generate HTML report
    GENERATE_HTML_REPORT (
        sorted_data,
        jinja_templates_list.first(),   // Pass the list of Jinja templates
        samplesheet.first(),            // Channel of one
        params_file.first(),            // Channel of one
        paramsJson,                     // JSON string can be used multiple times
        css_files_list.first()          // Channel of one (CSS files)
    )
    ch_versions = ch_versions.mix(GENERATE_HTML_REPORT.out.versions)

    emit:
    report    = GENERATE_HTML_REPORT.out.report
    versions  = ch_versions
}
