/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ASCC_GENOMIC      as GENOMIC      } from './ascc_genomic'
include { ASCC_ORGANELLAR   as ORGANELLAR   } from './ascc_organellar'

include { DECONTAMINATE_FASTA               } from '../subworkflows/local/decontaminate_fasta'

include { softwareVersionsToYAML            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText            } from '../subworkflows/local/utils_nfcore_ascc_pipeline'

// TODO: WHERE IS THIS FOR?
include { paramsSummaryMap                  } from 'plugin/nf-schema'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow ASCC {

    take:
    genomic_genomes         // channel: samplesheet read in from --input
    organellar_genomes      // channel: tuple(meta, reference)
    fcs_override            // bool
    genomic_fcs_samplesheet //
    organellar_fcs_samplesheet
    fcs_db                  // [path(path)]
    collected_reads         //
    scientific_name         // val(name)
    pacbio_database         // tuple [[meta.id], pacbio_database]
    ncbi_taxonomy_path
    ncbi_ranked_lineage_path
    nt_database_path
    diamond_nr_db_path
    diamond_uniprot_db_path
    taxid                   // val(taxid)
    nt_kraken_db_path
    vecscreen_database_path
    reads_path
    reads_layout
    reads_type
    btk_lineages
    btk_lineages_path
    barcodes

    main:
    ch_versions     = Channel.empty()

    //
    // WORKFLOW: Run main workflow for GENOMIC samples
    //
    //
    // TODO: FCS OVERRIDE VALUES NOW AVAILABLE
    GENOMIC (
        genomic_genomes,
        organellar_genomes,
        fcs_override,
        genomic_fcs_samplesheet,
        fcs_db,
        collected_reads,
        scientific_name,
        pacbio_database,
        ncbi_taxonomy_path,
        ncbi_ranked_lineage_path,
        nt_database_path,
        diamond_nr_db_path,
        diamond_uniprot_db_path,
        taxid,
        nt_kraken_db_path,
        vecscreen_database_path,
        reads_path,
        reads_layout,
        reads_type,
        btk_lineages,
        btk_lineages_path,
        barcodes
    )
    ch_versions         = ch_versions.mix(GENOMIC.out.versions)


    //
    // WORKFLOW: Run main workflow for ORGANELLAR samples
    //
    if ( !params.genomic_only ) {
        ORGANELLAR (
            organellar_genomes,
            fcs_override,
            organellar_fcs_samplesheet,
            fcs_db,
            collected_reads,
            scientific_name,
            pacbio_database,
            ncbi_taxonomy_path,
            ncbi_ranked_lineage_path,
            nt_database_path,
            diamond_nr_db_path,
            diamond_uniprot_db_path,
            taxid,
            nt_kraken_db_path,
            vecscreen_database_path,
            reads_path,
            reads_type,
            barcodes
        )
        ch_versions     = ch_versions.mix(ORGANELLAR.out.versions)
    }


    //
    // SUBWORKFLOW: 
    //



    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'ascc_software_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

workflow.onComplete {
    if (workflow.success) {
        try {
            def completionFile = file("${params.outdir}/workflow_completed.txt")
            completionFile.text = """
                Workflow completed successfully!
                Completed at: ${workflow.complete}
                Duration: ${workflow.duration}
                Success: ${workflow.success}
                Work directory: ${workflow.workDir}
                Exit status: ${workflow.exitStatus}
                Run name: ${workflow.runName}
                Session ID: ${workflow.sessionId}
                Project directory: ${workflow.projectDir}
                Launch directory: ${workflow.launchDir}
                Command line: ${workflow.commandLine}
            """.stripIndent()
            log.info "Completion file created: ${completionFile}"
        } catch (Exception e) {
            log.warn "Failed to create completion file: ${e.message}"
        }
    }
}
