/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ASCC GENOMIC REPORTING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CREATE_BTK_DATASET                            } from '../modules/local/blobtoolkit/create_dataset/main'
include { MERGE_BTK_DATASETS                            } from '../modules/local/blobtoolkit/merge_dataset/main'
include { ASCC_MERGE_TABLES                             } from '../modules/local/ascc/merge_tables/main'
include { AUTOFILTER_AND_CHECK_ASSEMBLY                 } from '../modules/local/autofilter/autofilter/main'
include { SANGER_TOL_BTK                                } from '../modules/local/sanger-tol/btk/main'
include { GENERATE_SAMPLESHEET                          } from '../modules/local/blobtoolkit/generate_samplesheet/main'
include { GENERATE_HTML_REPORT_WORKFLOW                 } from '../subworkflows/local/generate_html_report/main'

workflow ASCC_GENOMIC_REPORTING {
    take:
    reference_tuple_from_GG      // Channel: reference tuple from GENERATE_GENOME
    ej_dot_genome                // Channel: dot_genome from ESSENTIAL_JOBS
    ej_gc_coverage               // Channel: gc_content_txt from ESSENTIAL_JOBS
    include_workflow_steps       // List of included workflow steps
    exclude_workflow_steps       // List of excluded workflow steps
    ch_kmers                     // Channel: kmers from ASCC_GENOMIC_ANALYSIS
    ch_tiara                     // Channel: tiara from ASCC_GENOMIC_ANALYSIS
    ch_nt_blast                  // Channel: nt_blast from ASCC_GENOMIC_ANALYSIS
    ch_blast_lineage             // Channel: blast_lineage from ASCC_GENOMIC_ANALYSIS
    ch_btk_format                // Channel: btk_format from ASCC_GENOMIC_ANALYSIS
    nr_full                      // Channel: nr_full from ASCC_GENOMIC_ANALYSIS
    nr_hits                      // Channel: nr_hits from ASCC_GENOMIC_ANALYSIS
    un_full                      // Channel: un_full from ASCC_GENOMIC_ANALYSIS
    un_hits                      // Channel: un_hits from ASCC_GENOMIC_ANALYSIS
    ch_fcsgx                     // Channel: fcsgx from ASCC_GENOMIC_ANALYSIS
    ch_coverage                  // Channel: coverage from ASCC_GENOMIC_ANALYSIS
    ch_bam                       // Channel: bam from ASCC_GENOMIC_ANALYSIS
    ch_vecscreen                 // Channel: vecscreen from ASCC_GENOMIC_ANALYSIS
    ch_kraken1                   // Channel: kraken1 from ASCC_GENOMIC_ANALYSIS
    ch_kraken2                   // Channel: kraken2 from ASCC_GENOMIC_ANALYSIS
    ch_kraken3                   // Channel: kraken3 from ASCC_GENOMIC_ANALYSIS
    ch_fcsadapt                  // Channel: fcsadapt from ASCC_GENOMIC_ANALYSIS
    scientific_name              // val(name)
    btk_busco_run_mode           // val(mode)
    trailing_ns_report           // Channel: trailing_ns_report from ESSENTIAL_JOBS
    filter_fasta_sanitation_log  // Channel: filter_fasta_sanitation_log from ESSENTIAL_JOBS
    filter_fasta_length_filtering_log // Channel: filter_fasta_length_filtering_log from ESSENTIAL_JOBS
    pacbio_barcode_check_filtered // Channel: filtered from PACBIO_BARCODE_CHECK
    fcsgx_report                 // Channel: fcsgx_report from RUN_FCSGX
    fcsgx_taxonomy_report        // Channel: taxonomy_report from RUN_FCSGX
    ch_kmers_results             // Channel: kmers_results from ASCC_GENOMIC_ANALYSIS

    main:
    ch_versions = Channel.empty()
    
    // Initialize output channels
    create_summary = Channel.of([[],[]])
    busco_merge_btk = Channel.of([[],[]])
    ch_autofilt_assem = Channel.of([])
    ch_autofilt_indicator = Channel.of([])
    btk_bool = Channel.empty()

    //
    // BTK DATASET CREATION
    //
    if ( (include_workflow_steps.contains('create_btk_dataset') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("create_btk_dataset")
    ) {
        //
        // LOGIC: FOUND RACE CONDITION EFFECTING LONG RUNNING JOBS
        //          AND INPUT TO HERE ARE NOW MERGED AND MAPPED
        //          EMPTY CHANNELS ARE CHECKED AND DEFAULTED TO [[],[]]
        //
        ch_genomic_cbtk_input = reference_tuple_from_GG
            .map{ it -> tuple([
                id: it[0].id,
                taxid: it[0].taxid,
                sci_name: it[0].sci_name,
                process: "REFERENCE"], it[1])
            }
            .mix(
                ej_dot_genome.map{ it -> tuple([id: it[0].id, process: "GENOME"], it[1])},
                ch_kmers,
                ch_tiara,
                ch_nt_blast,
                // Use the BLAST top hits for BTK
                ch_btk_format,
                ch_fcsgx,
                ch_bam,
                ch_coverage,
                ch_kraken1,
                ch_kraken2,
                ch_kraken3,
                nr_full,
                un_full
            )
            .map { meta, file ->
                [meta.id, [meta: meta, file: file]]
            }
            .filter { id, data -> id != [] }
            .groupTuple()
            .map { id, data ->
                [id: id, data: data]
            }

        //
        // LOGIC: LIST OF PROCESSES TO CHECK FOR
        //
        def processes = [
            'REFERENCE', 'NT-BLAST', 'TIARA', 'Kraken 2', 'GENOME', 'KMERS',
            'FCSGX result', 'NR-FULL', 'UN-FULL', 'Mapped Bam', 'Coverage',
            'Kraken 1', 'Kraken 3'
        ]

        //
        // LOGIC: Create a channel for each process
        //
        def processChannels = processes.collectEntries { process ->
            [(process): ch_genomic_cbtk_input
                .map { sample ->
                    def data = sample.data.find { it.meta.process == process }
                    data ? [sample.id, data.meta, data.file] : [sample.id, [process: process], []]
                }
            ]
        }

        //
        // LOGIC: Combine all channels using a series of combine operations
        //
        def combined_channel = processChannels['REFERENCE']
        processes.tail().each { process ->
            combined_channel = combined_channel.combine(processChannels[process], by: 0)
        }

        //
        // MODULE: CREATE A BTK COMPATIBLE DATASET FOR NEW DATA
        //
        CREATE_BTK_DATASET (
            combined_channel,
            Channel.fromPath(params.ncbi_taxonomy_path).first(),
            scientific_name
        )
        ch_versions             = ch_versions.mix(CREATE_BTK_DATASET.out.versions)

        create_summary          = CREATE_BTK_DATASET.out.create_summary.map{ it -> tuple([id: it[0].id, process: "C_BTK_SUM"], it[1])}
    }

    //
    // AUTOFILTER ASSEMBLY
    //
    if (
        ( include_workflow_steps.contains('tiara') && include_workflow_steps.contains('fcs-gx') && include_workflow_steps.contains("autofilter_assembly") && !exclude_workflow_steps.contains("autofilter_assembly") ) ||
        ( include_workflow_steps.contains('ALL') && !exclude_workflow_steps.contains("autofilter_assembly") )
    ) {
        //
        // LOGIC: FILTER THE INPUT FOR THE AUTOFILTER STEP
        //          - We can't just combine on meta.id as some of the Channels have other data
        //              in there too so we just sanitise, and _then_ combine on 0, and
        //              _then_ add back in the taxid as we need that for this process.
        //              Thankfully taxid is a param so easy enough to add back in.
        //                  Actually, it just makes more sense to passs in as its own channel.
        //

        autofilter_input_formatted = reference_tuple_from_GG
            .map{ it -> tuple([id: it[0].id], it[1])}
            .combine(
                ch_tiara
                    .map{ it -> tuple([id: it[0].id], it[1])},
                by: 0
            )
            .combine(
                ch_fcsgx
                    .map{ it -> tuple([id: it[0].id], it[1])},
                by: 0
            )
            .combine(
                Channel.fromPath(params.ncbi_ranked_lineage_path)
            )
            .combine(
                Channel.of(params.taxid)
            )
            .multiMap{
                meta, ref, tiara, fcs, ncbi, thetaxid ->
                    reference:  tuple([id: meta.id, taxid: thetaxid], ref)
                    tiara_file: tuple(meta, tiara)
                    fcs_file:   tuple(meta, fcs)
                    ncbi_rank:  ncbi
            }

        //
        // MODULE: AUTOFILTER ASSEMBLY BY TIARA AND FCSGX RESULTS
        //
        AUTOFILTER_AND_CHECK_ASSEMBLY (
            autofilter_input_formatted.reference,
            autofilter_input_formatted.tiara_file,
            autofilter_input_formatted.fcs_file,
            autofilter_input_formatted.ncbi_rank
        )
        ch_autofilt_assem       = AUTOFILTER_AND_CHECK_ASSEMBLY.out.decontaminated_assembly.map{it[1]}
        ch_autofilt_indicator   = AUTOFILTER_AND_CHECK_ASSEMBLY.out.indicator_file

        //
        // LOGIC: BRANCH THE CHANNEL ON WHETHER OR NOT THERE IS ABNORMAL CONTAMINATION IN THE
        //          OUTPUT FILE.
        //
        btk_bool = AUTOFILTER_AND_CHECK_ASSEMBLY.out.alarm_file
            .map { file -> file.text.trim() }
            .branch { it ->
                run_btk     : it.contains("YES_ABNORMAL_CONTAMINATION")
                dont_run    : true
            }

        ch_versions             = ch_versions.mix(AUTOFILTER_AND_CHECK_ASSEMBLY.out.versions)
    }

    //
    // BTK/BUSCO PROCESSING
    //
    if (
        (
            !exclude_workflow_steps.contains("btk_busco") &&
            ((include_workflow_steps.contains('btk_busco') && include_workflow_steps.contains("autofilter_assembly")) || include_workflow_steps.contains('ALL')) &&
            btk_busco_run_mode == "conditional" &&
            btk_bool.run_btk
        ) ||
        (
            !exclude_workflow_steps.contains("btk_busco") &&
            ((include_workflow_steps.contains('btk_busco') && include_workflow_steps.contains("autofilter_assembly") && include_workflow_steps.contains("create_btk_dataset")) || include_workflow_steps.contains('ALL')) &&
            btk_busco_run_mode == "mandatory"
        )
    ) {
        //
        // MODULE: THIS MODULE FORMATS THE INPUT DATA IN A SPECIFIC CSV FORMAT FOR
        //          USE IN THE BTK PIPELINE
        //
        GENERATE_SAMPLESHEET (
            reference_tuple_from_GG,
            params.reads_path,
            Channel.of(params.reads_layout),
            AUTOFILTER_AND_CHECK_ASSEMBLY.out.alarm_file
        )
        ch_versions         = ch_versions.mix(GENERATE_SAMPLESHEET.out.versions)

        //
        // LOGIC: STRIP THE META DATA DOWN TO id AND COMBINE ON THAT.
        //
        coverage_id = GENERATE_SAMPLESHEET.out.csv
            .map{ meta, csv ->
                tuple(
                    [ id: meta.id ],
                    csv
                )
            }

        combined_input = reference_tuple_from_GG
            .map{ it -> tuple([id:it[0].id], it[1])}
            .combine(coverage_id, by: 0)
            .multiMap { meta_1, ref, csv ->
                reference: [meta_1, ref]
                samplesheet: csv
            }

        //
        // PIPELINE: PREPARE THE DATA FOR USE IN THE SANGER-TOL/BLOBTOOLKIT PIPELINE
        //              WE ARE USING THE PIPELINE HERE AS A MODULE THIS REQUIRES IT
        //              TO BE USED AS A AN INTERACTIVE JOB ON WHAT EVER EXECUTOR YOU ARE USING.
        //              This will also eventually check for the above run_btk boolean from
        //              autofilter
        SANGER_TOL_BTK (
            combined_input.reference,
            combined_input.samplesheet,
            params.diamond_uniprot_database_path,
            params.nt_database_path,
            params.diamond_uniprot_database_path,
            params.ncbi_taxonomy_path,
            params.reads_path,
            params.busco_lineages_folder,
            params.busco_lineages,
            params.taxid,
        )
        ch_versions             = ch_versions.mix(SANGER_TOL_BTK.out.versions)

        //
        // MODULE: MERGE THE TWO BTK FORMATTED DATASETS INTO ONE DATASET FOR EASIER USE
        //
        merged_channel = CREATE_BTK_DATASET.out.btk_datasets
            .map { meta, file -> [meta.id, [meta, file]] }
            .join(
                SANGER_TOL_BTK.out.dataset
                    .map { meta, file ->
                        [meta.id, [meta, file]]
                })
            .map { id, ref, btk -> [ref[0], ref[1], btk[1]] }

        MERGE_BTK_DATASETS (
            merged_channel
        )
        ch_versions             = ch_versions.mix(MERGE_BTK_DATASETS.out.versions)
        busco_merge_btk         = MERGE_BTK_DATASETS.out.busco_summary_tsv
    }

    //
    // TABLE MERGING
    //
    if (
        !exclude_workflow_steps.contains("essentials") && !exclude_workflow_steps.contains("merge") && !exclude_workflow_steps.contains("ALL")
    ) {
        //
        // LOGIC: FOUND RACE CONDITION EFFECTING LONG RUNNING JOBS
        //          AND INPUT TO HERE ARE NOW MERGED AND MAPPED
        //          EMPTY CHANNELS ARE CHECKED AND DEFAULTED TO [[],[]]
        //
        ascc_merged_data = ej_gc_coverage
            .map{ it -> tuple([
                id: it[0].id,
                process: "GC_COV"], it[1])
            }
            .mix(
                ej_dot_genome.map{ it -> tuple([id: it[0].id, process: "GENOME"], it[1])},
                create_summary,
                busco_merge_btk.map{ it -> tuple([id: it[0].id, process: "BUSCO_MERGE"], it[1])},
                ch_kmers,
                ch_tiara,
                ch_fcsgx,
                ch_coverage,
                ch_kraken3,
                ch_blast_lineage,
                nr_hits,
                un_hits
            )
            .map { meta, file ->
                [meta.id, [meta: meta, file: file]]
            }
            .filter { id, data -> id != [] }
            .groupTuple()
            .map { id, data ->
                [id: id, data: data]
            }

        def processes = [
            'GC_COV', 'Coverage', 'TIARA',
            'Kraken 3', 'NT-BLAST-LINEAGE', 'KMERS', 'NR-HITS', 'UN-HITS',
            'C_BTK_SUM', 'BUSCO_MERGE','FCSGX result'
        ]

        def processChannels = processes.collectEntries { process ->
            [(process): ascc_merged_data
                .map { sample ->
                    def data = sample.data.find { it.meta.process == process }
                    data ? [sample.id, data.meta, data.file] : [sample.id, [process: process], []]
                }
            ]
        }

        def ascc_combined_channels = processChannels['GC_COV']
        processes.tail().each { process ->
            ascc_combined_channels = ascc_combined_channels
                                    .combine(processChannels[process], by: 0)
        }

        //
        // SUBWORKFLOW: MERGES DATA THAT IS NOT USED IN THE CREATION OF THE BTK_DATASETS FOLDER
        //
        ASCC_MERGE_TABLES (
            ascc_combined_channels.map { it[1..-1] } // Remove the first item in tuple (mapping key)
        )
        ch_versions             = ch_versions.mix(ASCC_MERGE_TABLES.out.versions)
    }

    //
    // HTML REPORT GENERATION
    //
    if ( (include_workflow_steps.contains('html_report') || include_workflow_steps.contains('ALL')) &&
            !exclude_workflow_steps.contains("html_report")
    ) {
        // Create channels for HTML report inputs
        // Initialize all channels as empty by default
        local_empty_barcode_channel = Channel.of([[id: "empty"],[]])  // Renamed for clarity
        ch_fcs_adaptor_euk = Channel.of([[id: "empty"],[]])
        ch_fcs_adaptor_prok = Channel.of([[id: "empty"],[]])
        ch_trim_ns_results = Channel.of([[id: "empty"],[]])
        ch_vecscreen_results = Channel.of([[id: "empty"],[]])
        local_empty_autofilter_channel = Channel.of([[id: "empty"],[]])  // Renamed for clarity
        ch_merged_table = Channel.of([[id: "empty"],[]])
        ch_phylum_counts = Channel.of([[id: "empty"],[]])  // Initialize phylum counts channel
        local_empty_kmers_channel = Channel.of([[id: "empty"],[]]) // Renamed for clarity
        ch_fasta_sanitation_log = Channel.of([[id: "empty"],[]])
        ch_fasta_length_filtering_log = Channel.of([[id: "empty"],[]])
        ch_fcsgx_report_txt = Channel.of([[id: "empty"],[]])
        ch_fcsgx_taxonomy_rpt = Channel.of([[id: "empty"],[]])
        
        // Determine which barcode channel to pass based on the condition
        def final_barcode_channel_for_report
        if ((include_workflow_steps.contains('pacbio_barcodes') || include_workflow_steps.contains('ALL')) &&
                !exclude_workflow_steps.contains("pacbio_barcodes")) {
            // If pacbio_barcodes step included, use the channel from the input
            final_barcode_channel_for_report = pacbio_barcode_check_filtered
        } else {
            // Otherwise, use the initialized empty channel
            final_barcode_channel_for_report = local_empty_barcode_channel
        }
        
        if ((include_workflow_steps.contains('fcs-adaptor') || include_workflow_steps.contains('ALL')) &&
                !exclude_workflow_steps.contains("fcs-adaptor")) {
            ch_fcs_adaptor_euk = ch_fcsadapt.map { meta, file1, file2 -> [meta, file1] }
            ch_fcs_adaptor_prok = ch_fcsadapt.map { meta, file1, file2 -> [meta, file2] }
        }
        
        ch_trim_ns_results = trailing_ns_report
        
        if ((include_workflow_steps.contains('vecscreen') || include_workflow_steps.contains('ALL')) &&
                !exclude_workflow_steps.contains("vecscreen")) {
            ch_vecscreen_results = ch_vecscreen
        }
        
        // Determine which autofilter channel to pass based on the condition
        def final_autofilter_channel_for_report
        if ((include_workflow_steps.contains('autofilter_assembly') || include_workflow_steps.contains('ALL')) &&
                !exclude_workflow_steps.contains("autofilter_assembly")) {
            // If autofilter step included, use the fcs_tiara_summary file (CSV with contamination results)
            // This contains the actual autofiltering results, not just an indicator file
            final_autofilter_channel_for_report = AUTOFILTER_AND_CHECK_ASSEMBLY.out.fcs_tiara_summary
        } else {
            // Otherwise, use the initialized empty channel
            final_autofilter_channel_for_report = local_empty_autofilter_channel
        }
        
        // Get the kmers results if the kmers workflow was run
        if ((include_workflow_steps.contains('kmers') || include_workflow_steps.contains('ALL')) &&
                !exclude_workflow_steps.contains("kmers")) {
            // Use the kmers_results channel that contains the results directories with PNG files and metrics
            // The channel passed in via 'take' will be used if this condition is true
        }

        // Determine which kmer channel to pass based on the condition
        def final_kmers_channel_for_report
        if ((include_workflow_steps.contains('kmers') || include_workflow_steps.contains('ALL')) &&
                !exclude_workflow_steps.contains("kmers")) {
            // If kmers step included, use the channel passed via the 'take:' block
            final_kmers_channel_for_report = ch_kmers_results // This refers to the input channel from 'take:'
        } else {
            // Otherwise, use the initialized empty channel
            final_kmers_channel_for_report = local_empty_kmers_channel
        }
        
        // Get the merged table and phylum counts if the merge workflow was run
        if (!exclude_workflow_steps.contains("essentials") && !exclude_workflow_steps.contains("merge")) {
            ch_merged_table = ASCC_MERGE_TABLES.out.merged_table
            ch_phylum_counts = ASCC_MERGE_TABLES.out.phylum_counts
        }
        
        // Get the FASTA sanitation log if the filter_fasta process was run
        ch_fasta_sanitation_log = filter_fasta_sanitation_log
        ch_fasta_length_filtering_log = filter_fasta_length_filtering_log
        
        // Capture FCS-GX report files if available
        if ((include_workflow_steps.contains('fcs-gx') || include_workflow_steps.contains('ALL')) &&
                !exclude_workflow_steps.contains("fcs-gx")) {
            ch_fcsgx_report_txt = fcsgx_report
            ch_fcsgx_taxonomy_rpt = fcsgx_taxonomy_report
        }
        
        // Create channels for the input samplesheet and YAML parameters file
        ch_samplesheet_path = Channel.fromPath(params.input)
        
        // Handle the params_file parameter - check if it exists and create a channel
        if (params.containsKey('params_file') && params.params_file) {
            ch_params_file = Channel.fromPath(params.params_file)
        } else {
            ch_params_file = Channel.value([])
        }
        
        // Add debug logging
        log.info "HTML Report Generation - Input Channels in genomic workflow:"
        log.info "final_barcode_channel_for_report: ${final_barcode_channel_for_report.dump()}"
        log.info "ch_fcs_adaptor_euk: ${ch_fcs_adaptor_euk.dump()}"
        log.info "ch_fcs_adaptor_prok: ${ch_fcs_adaptor_prok.dump()}"
        log.info "ch_trim_ns_results: ${ch_trim_ns_results.dump()}"
        log.info "ch_vecscreen_results: ${ch_vecscreen_results.dump()}"
        log.info "final_autofilter_channel_for_report: ${final_autofilter_channel_for_report.dump()}"
        log.info "ch_fasta_sanitation_log: ${ch_fasta_sanitation_log.dump()}"
        log.info "ch_fasta_length_filtering_log: ${ch_fasta_length_filtering_log.dump()}"
        log.info "ch_phylum_counts: ${ch_phylum_counts.dump()}"
        log.info "ch_samplesheet_path: ${ch_samplesheet_path.dump()}"
        log.info "ch_params_file: ${ch_params_file.dump()}"
        
        // Get all Jinja templates and CSS files
        ch_jinja_templates_list = Channel.fromPath("${baseDir}/assets/templates/*.jinja").collect()
        ch_css_files = Channel.fromPath("${baseDir}/assets/css/*.css").collect()

        // Create a channel with the reference file
        ch_reference_file = reference_tuple_from_GG

        // Determine which BTK dataset channel to pass based on the condition
        def final_btk_dataset_channel_for_report
        if ((include_workflow_steps.contains('create_btk_dataset') || include_workflow_steps.contains('ALL')) &&
                !exclude_workflow_steps.contains("create_btk_dataset")) {
            // If create_btk_dataset step included, use the channel from CREATE_BTK_DATASET
            final_btk_dataset_channel_for_report = CREATE_BTK_DATASET.out.btk_datasets
        } else {
            // Otherwise, use an empty channel
            final_btk_dataset_channel_for_report = Channel.of([[id: "empty"],[]])
        }

        GENERATE_HTML_REPORT_WORKFLOW (
            final_barcode_channel_for_report,  // Pass the determined channel here
            ch_fcs_adaptor_euk,
            ch_fcs_adaptor_prok,
            ch_trim_ns_results,
            ch_vecscreen_results,
            final_autofilter_channel_for_report,  // Pass the determined channel here
            ch_merged_table,
            ch_phylum_counts,         // Pass the phylum counts channel
            final_kmers_channel_for_report, // Pass the determined channel here
            ch_reference_file,
            ch_fasta_sanitation_log,
            ch_fasta_length_filtering_log,
            ch_jinja_templates_list, // Pass the list of Jinja templates
            ch_samplesheet_path,
            ch_params_file,
            ch_fcsgx_report_txt,      // Pass FCS-GX report txt
            ch_fcsgx_taxonomy_rpt,    // Pass FCS-GX taxonomy rpt
            final_btk_dataset_channel_for_report, // Pass BTK dataset
            ch_css_files              // Pass CSS files
        )
        ch_versions = ch_versions.mix(GENERATE_HTML_REPORT_WORKFLOW.out.versions)
    }

    emit:
    versions = ch_versions
}
