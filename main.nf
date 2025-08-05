#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    sanger-tol/ascc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/sanger-tol/ascc
    Website: https://pipelines.tol.sanger.ac.uk/ascc
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ASCC                      } from './workflows/ascc'

include { PIPELINE_INITIALISATION   } from './subworkflows/local/utils_nfcore_ascc_pipeline'
include { PIPELINE_COMPLETION       } from './subworkflows/local/utils_nfcore_ascc_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow SANGERTOL_ASCC {

    take:
    genomic                     // Genomic fasta tuples
    organelles                  // Organellar fasta tuples
    fcs_override                // bool
    genomic_fcs_samplesheet     // genomic data tuple from fcs_override samplesheet generated from external fcs process
    organellar_fcs_samplesheet  // organellar data tuple from fcs_override samplesheet generated from external fcs process
    fcs_db                      // [path(path)]
    read_files                  // Read files
    scientific_name             // Scientific name
    pacbio_db                   // Pacbio database
    ncbi_taxonomy_path          // NCBI taxonomy path
    ncbi_ranked_lineage_path    // ranked_lineage file path
    nt_database_path            // NT database path
    diamond_nr_db_path          // NR database path
    diamond_uniprot_db_path     // UNIPROT database path
    taxid                       // taxid of input assembly
    nt_kraken_db_path           // path to custom kraken database
    vecscreen_database_path     // Path to custom Vecscreen database
    reads_path                  // tuple of paths to input read files
    reads_layout                // Layout of input reads [SINGLE, PAIRED]
    reads_type                  // Type of input reads [HIFI, CLR, ILLUMINA, etc...]
    btk_lineages                // [1_odb10, 2_odb10] ODB10 lineages to run against the assembly
    btk_lineages_path           // Path to the ODB10 lineages

    main:

    //
    // WORKFLOW: Run pipeline
    //
    ASCC (
        genomic,
        organelles,
        fcs_override,
        genomic_fcs_samplesheet,
        organellar_fcs_samplesheet,
        fcs_db,
        read_files,
        scientific_name,
        pacbio_db,
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
        btk_lineages_path
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )


    //
    // WORKFLOW: MAIN ASCC WORKFLOW FILE THAT SEPERATES INTO GENOMIC AND ORGANELLAR
    //
    SANGERTOL_ASCC (
        PIPELINE_INITIALISATION.out.main_genomes,
        PIPELINE_INITIALISATION.out.organellar_genomes,
        params.fcs_override,
        PIPELINE_INITIALISATION.out.geno_fcs_samplesheet,
        PIPELINE_INITIALISATION.out.orga_fcs_samplesheet,
        PIPELINE_INITIALISATION.out.fcs_gx_database,
        PIPELINE_INITIALISATION.out.collected_reads,
        params.scientific_name,
        PIPELINE_INITIALISATION.out.pacbio_db,
        Channel.fromPath(params.ncbi_taxonomy_path),
        Channel.fromPath(params.ncbi_ranked_lineage_path),
        Channel.fromPath(params.nt_database_path),
        Channel.fromPath(params.diamond_nr_database_path),
        Channel.fromPath(params.diamond_uniprot_database_path),
        Channel.of(params.taxid),
        Channel.fromPath(params.nt_kraken_database_path),
        Channel.fromPath(params.vecscreen_database_path),
        Channel.from(params.reads_path),
        Channel.of(params.reads_layout),
        Channel.of(params.reads_type),
        Channel.of(params.busco_lineages),
        Channel.fromPath(params.busco_lineages_folder)
    )


    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        []
    )
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
