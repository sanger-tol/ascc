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
include { NEXTFLOW_RUN as SANGER_TOL_BTK_CASCADE        } from '../modules/local/run/main'
include { GENERATE_HTML_REPORT_WORKFLOW                 } from '../subworkflows/local/generate_html_report/main'

workflow ASCC_GENOMIC_REPORTING {
    take:
    reference_tuple_from_GG      // Channel: reference tuple from GENERATE_GENOME
    ej_dot_genome                // Channel: dot_genome from ESSENTIAL_JOBS
    ej_gc_coverage               // Channel: gc_content_txt from ESSENTIAL_JOBS
    kmers                        // Channel: from analysis workflow
    tiara                        // Channel: from analysis workflow
    nt_blast                     // Channel: from analysis workflow
    blast_lineage                // Channel: from analysis workflow
    btk_format                   // Channel: from analysis workflow
    nr_full                      // Channel: from analysis workflow
    nr_hits                      // Channel: from analysis workflow
    un_full                      // Channel: from analysis workflow
    un_hits                      // Channel: from analysis workflow
    fcsgx                        // Channel: from analysis workflow
    coverage                     // Channel: from analysis workflow
    bam                          // Channel: from analysis workflow
    vecscreen                    // Channel: from analysis workflow
    kraken1                      // Channel: from analysis workflow
    kraken2                      // Channel: from analysis workflow
    kraken3                      // Channel: from analysis workflow
    fcsadapt                     // Channel: from analysis workflow
    scientific_name              // val(name)
    btk_busco_run_mode          // val(mode)
    trailing_ns_report          // Channel: from essential jobs
    filter_fasta_sanitation_log // Channel: from essential jobs
    filter_fasta_length_filtering_log // Channel: from essential jobs
    pacbio_barcode_check_filtered // Channel: from analysis workflow
    fcsgx_report                // Channel: from analysis workflow
    fcsgx_taxonomy_report       // Channel: from analysis workflow
    kmers_results               // Channel: from analysis workflow

    main:
    ch_versions = Channel.empty()

    // Initialize output channels
    create_summary = Channel.empty()
    busco_merge_btk = Channel.empty()
    ch_autofilt_assem = Channel.empty()
    ch_autofilt_indicator = Channel.empty()
    ch_phylum_counts = Channel.empty()

    //
    // BTK DATASET CREATION
    //
    if (params.run_create_btk_dataset == "both" || params.run_create_btk_dataset == "genomic") {
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
                kmers,
                tiara,
                nt_blast,
                btk_format,
                fcsgx,
                bam,
                coverage,
                kraken1,
                kraken2,
                kraken3,
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

        def processes = [
            'REFERENCE', 'NT-BLAST', 'TIARA', 'Kraken 2', 'GENOME', 'KMERS',
            'FCSGX result', 'NR-FULL', 'UN-FULL', 'Mapped Bam', 'Coverage',
            'Kraken 1', 'Kraken 3'
        ]

        def processChannels = processes.collectEntries { process ->
            [(process): ch_genomic_cbtk_input
                .map { sample ->
                    def data = sample.data.find { it.meta.process == process }
                    data ? [sample.id, data.meta, data.file] : [sample.id, [process: process], []]
                }
            ]
        }

        def combined_channel = processChannels['REFERENCE']
        processes.tail().each { process ->
            combined_channel = combined_channel.combine(processChannels[process], by: 0)
        }

        CREATE_BTK_DATASET (
            combined_channel,
            params.ncbi_taxonomy_path,
            scientific_name
        )
        ch_versions             = ch_versions.mix(CREATE_BTK_DATASET.out.versions)

        create_summary          = CREATE_BTK_DATASET.out.create_summary.map{ it -> tuple([id: it[0].id, process: "C_BTK_SUM"], it[1])}
    }

    //
    // AUTOFILTER ASSEMBLY
    //
    if ((params.run_tiara == "both" || params.run_tiara == "genomic") &&
        (params.run_fcsgx == "both" || params.run_fcsgx == "genomic") &&
        (params.run_autofilter_assembly == "both" || params.run_autofilter_assembly == "genomic")) {
        autofilter_input_formatted = reference_tuple_from_GG
            .map{ it -> tuple([id: it[0].id], it[1])}
            .combine(
                tiara
                    .map{ it -> tuple([id: it[0].id], it[1])},
                by: 0
            )
            .combine(
                fcsgx
                    .map{ it -> tuple([id: it[0].id], it[1])},
                by: 0
            )
            .combine(
                Channel.of(params.ncbi_ranked_lineage_path)
            )
            .combine(
                Channel.of(params.taxid)
            )
            .multiMap{
                meta, ref, tiara_file, fcs_file, ncbi, thetaxid ->
                    reference:  tuple([id: meta.id, taxid: thetaxid], ref)
                    tiara_file: tuple(meta, tiara_file)
                    fcs_file:   tuple(meta, fcs_file)
                    ncbi_rank:  ncbi
            }

        AUTOFILTER_AND_CHECK_ASSEMBLY (
            autofilter_input_formatted.reference,
            autofilter_input_formatted.tiara_file,
            autofilter_input_formatted.fcs_file,
            autofilter_input_formatted.ncbi_rank
        )
        ch_autofilt_assem       = AUTOFILTER_AND_CHECK_ASSEMBLY.out.decontaminated_assembly.map{it[1]}
        ch_autofilt_indicator   = AUTOFILTER_AND_CHECK_ASSEMBLY.out.indicator_file

        btk_bool = AUTOFILTER_AND_CHECK_ASSEMBLY.out.alarm_file
            .map { file -> file.text.trim() }
            .branch { it ->
                run_btk     : it.contains("YES_ABNORMAL_CONTAMINATION")
                dont_run    : true
            }

        ch_versions             = ch_versions.mix(AUTOFILTER_AND_CHECK_ASSEMBLY.out.versions)
    } else {
        ch_autofilt_assem       = Channel.empty()
        ch_autofilt_indicator   = Channel.empty()
        btk_bool = Channel.empty()
    }

    //
    // TABLE MERGING
    //
    if ((params.run_essentials == "both" || params.run_essentials == "genomic") &&
        (params.run_merge_datasets == "both" || params.run_merge_datasets == "genomic")) {
        ascc_merged_data = ej_gc_coverage
            .map{ it -> tuple([
                id: it[0].id,
                process: "GC_COV"], it[1])
            }
            .mix(
                ej_dot_genome.map{ it -> tuple([id: it[0].id, process: "GENOME"], it[1])},
                create_summary,
                busco_merge_btk.map{ it -> tuple([id: it[0].id, process: "BUSCO_MERGE"], it[1])},
                kmers,
                tiara,
                fcsgx,
                coverage,
                kraken3,
                blast_lineage,
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

        ASCC_MERGE_TABLES (
            ascc_combined_channels.map { it[1..-1] }
        )
        ch_versions             = ch_versions.mix(ASCC_MERGE_TABLES.out.versions)
        ch_phylum_counts = ASCC_MERGE_TABLES.out.phylum_counts
    }

    //
    // HTML REPORT GENERATION
    //
    if (params.run_html_report == "both" || params.run_html_report == "genomic") {
        
        // Get all Jinja templates and CSS files
        ch_jinja_templates_list = Channel.fromPath("${baseDir}/assets/templates/*.jinja").collect()
        ch_css_files = Channel.fromPath("${baseDir}/assets/css/*.css").collect()
        
        // Create samplesheet and params file channels
        samplesheet = Channel.fromPath(params.input)
        params_file = Channel.value(file("${workflow.launchDir}/nextflow.config"))
        
        // Create empty channel with proper meta structure for missing data
        empty_meta_channel = Channel.of([[id: "empty"], []])
        
        // Get BTK dataset if available
        btk_dataset_channel = (params.run_create_btk_dataset == "both" || 
                              params.run_create_btk_dataset == "genomic") ? 
                             CREATE_BTK_DATASET.out.btk_datasets : empty_meta_channel
        
        // Get merged table if available
        merged_table_channel = ((params.run_essentials == "both" || params.run_essentials == "genomic") &&
                               (params.run_merge_datasets == "both" || params.run_merge_datasets == "genomic")) ? 
                               ASCC_MERGE_TABLES.out.merged_table : empty_meta_channel
        
        // Get autofilter results if available
        autofilter_channel = ((params.run_tiara == "both" || params.run_tiara == "genomic") &&
                             (params.run_fcsgx == "both" || params.run_fcsgx == "genomic") &&
                             (params.run_autofilter_assembly == "both" || params.run_autofilter_assembly == "genomic")) ? 
                             AUTOFILTER_AND_CHECK_ASSEMBLY.out.fcs_tiara_summary : empty_meta_channel
        
        // Handle FCS adaptor channels properly
        if ((params.run_fcs_adaptor == "both" || params.run_fcs_adaptor == "genomic")) {
            fcs_adaptor_euk_channel = fcsadapt.map { meta, euk, prok -> [meta, euk] }
            fcs_adaptor_prok_channel = fcsadapt.map { meta, euk, prok -> [meta, prok] }
        } else {
            fcs_adaptor_euk_channel = empty_meta_channel
            fcs_adaptor_prok_channel = empty_meta_channel
        }
        
        GENERATE_HTML_REPORT_WORKFLOW (
            pacbio_barcode_check_filtered.ifEmpty(empty_meta_channel),  // barcode_results
            fcs_adaptor_euk_channel,        // fcs_adaptor_euk
            fcs_adaptor_prok_channel,       // fcs_adaptor_prok
            trailing_ns_report.ifEmpty(empty_meta_channel),             // trim_ns_results
            vecscreen.ifEmpty(empty_meta_channel),                      // vecscreen_results
            autofilter_channel,             // autofilter_results
            merged_table_channel,           // merged_table
            ch_phylum_counts.ifEmpty(empty_meta_channel),               // phylum_counts
            kmers_results.ifEmpty(empty_meta_channel),                  // kmers_results
            reference_tuple_from_GG,        // reference_fasta
            filter_fasta_sanitation_log.ifEmpty(empty_meta_channel),    // fasta_sanitation_log
            filter_fasta_length_filtering_log.ifEmpty(empty_meta_channel), // fasta_length_filtering_log
            ch_jinja_templates_list,        // jinja_templates_list
            samplesheet,                    // samplesheet
            params_file,                    // params_file
            fcsgx_report.ifEmpty(empty_meta_channel),                   // fcsgx_report_txt
            fcsgx_taxonomy_report.ifEmpty(empty_meta_channel),          // fcsgx_taxonomy_rpt
            btk_dataset_channel,            // btk_dataset
            ch_css_files                    // css_files_list
        )
        ch_versions = ch_versions.mix(GENERATE_HTML_REPORT_WORKFLOW.out.versions)
    }

    emit:
    versions = ch_versions
}
