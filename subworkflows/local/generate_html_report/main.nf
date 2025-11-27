include { GENERATE_HTML_REPORT } from '../../../modules/local/generate_html_report/main'

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
    ch_versions = channel.empty()

    // Convert params to JSON for passing to the HTML report
    def paramsJson = groovy.json.JsonOutput.toJson(params)

    //
    // LOGIC: COMBINE ALL INPUT CHANNELS AND ANNOTATE THEM WITH PROCESS TAGS
    //
    fcs_adaptor
        .multiMap { meta, file1, file2 ->
            euk:    [meta + [process: "FCS_ADAPTOR_EUK"],   file1]
            prok:   [meta + [process: "FCS_ADAPTOR_PROK"],  file2]
        }
        .set { fcs_adaptor_split }

    ch_fcs_euk  = fcs_adaptor_split.euk
                    .ifEmpty{ [[process: "FCS_ADAPTOR_EUK"], []] }

    ch_fcs_prok = fcs_adaptor_split.prok
                    .ifEmpty{ [[process: "FCS_ADAPTOR_PROK"], []] }


    barcode_results
        .mix(
            trim_ns_results,
            autofilter_results,
            merged_table,
            phylum_counts,
            ch_fcs_euk,
            ch_fcs_prok,
            vecscreen_results,
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
        // BELOW SHOULD FILTER AWAY [] and null VALUE CHANNELS
        // THIS GETS MADE NOW THAT ALL CHANNELS HAVE PROCESS TAGS
        // IT FORCES THE CREATION OF THE NULL CHANNEL
        .filter { id, data -> id }
        // GROUP ON 'KEY' AND RE-MAP INTO A MORE SENSIBLE STRUCTURE
        .groupTuple()
        .map { id, data ->
            [id: id, data: data]
        }
        .set { all_data }

        // LIST OF EXPECTED INPUT DATA
        def processes = [
            "REFERENCE",
            "BARCODES", "REFERENCE_SANI_LOG", "REFERENCE_FILT_LOG",
            "TRAILING_NS", "FCSGX_REPORT", "FCSGX_TAX_REPORT",
            "VECSCREEN", "KMER_RESULTS", "FCS_ADAPTOR_EUK",
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
        def combined_channels = processChannels["REFERENCE"]

        processes.tail().each { process ->
            combined_channels = combined_channels
                .combine(processChannels[process], by: 0)
        }

        // ABOVE PROCESSES ARE NOT DETEMINISTIC, HOWEVER, WE SHOULD BE ABLE TO
        // SORT THEM INTO THE GIVEN ORDER
        combined_channels
            .map { sample ->
                def id = sample[0]  // CREATE COMMON KEY
                def dataItems = sample[1..-1]  // DATA CHANNELS

                // CREATE A MAP ON THE PROCESSES LIST, THIS WILL BE TRUTH ORDER
                def processOrder = processes.withIndex().collectEntries { proc, idx ->
                    [proc, idx]
                }

                // SORT DATA [meta, file] BY THE PROCESS ORDER
                //
                def metaFilePairs = dataItems.collate(2)
                def sortedData =  metaFilePairs
                    .sort { metaFile ->
                        processOrder[metaFile[0].process]
                    }
                    .collect()

                // CHANNEL OUTPUTS AS [meta,[meta, file]]
                [[ id: id ], sortedData]
            }
            .map { item ->
                // FLATTEN INTO [meta0, meta1, file1, meta2, file2...] KEEPING []
                // THIS ALLOWS US TO USE INDEXES LATER ON
                //
                // START WITH SEED META
                def result = [item[0]]

                item[1].each { pair ->
                    // SUB-META
                    result.add(pair[0])

                    // FILE or [], [] MUST BE PRESERVED
                    result.add(pair[1])
                }

                return result
            }
            // SO WE HAVE THE SAME NUMBER OF INPUT CHANNELS TO PROCESS
            .combine ( samplesheet )
            .combine ( params_file.map { [it] } )
            .combine ( channel.value(paramsJson) )

            // CAN'T NAME THE ITEMS, SIMPLY TOO MANY TO INDEX IT IS
            .multiMap { item ->
                // META TO ACT AS SEED FOR LIST EXTENSION
                def result = [item[0]]

                // DATA FILES END UP NESTED, THIS PULLS THEM OUT OF THAT.
                // USING -INDEX SHOULD MEAN IT'S OK TO ADD MORE IN FUTURE
                // AS LONG AS IT IS IN THE PROCESS MAP... AND THE MODULE
                def data_files = item[1..-4]
                    .each { iter ->
                        result.add(iter)
                    }

                data: result
                samplesheet: item[-3]
                params: item[-2]
                json: item[-1]
            }
            .set { sorted_data }


    //
    // MODULE: GENERATE A HTML REPORT FOR END USER
    //
    GENERATE_HTML_REPORT (
        sorted_data.data,
        channel.fromPath("${projectDir}/assets/templates/*.jinja").collect(),   // Pass the list of Jinja templates
        sorted_data.samplesheet,                                                // channel.of one
        sorted_data.params,                                                     // channel.of one
        sorted_data.json,                                                       // JSON string can be used multiple times
        channel.fromPath("${projectDir}/assets/css/*.css").collect()            // channel.of one (CSS files)
    )
    ch_versions = ch_versions.mix(GENERATE_HTML_REPORT.out.versions)

    emit:
    report    = GENERATE_HTML_REPORT.out.report
    versions  = ch_versions
}
