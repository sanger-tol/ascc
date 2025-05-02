/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ASCC_GENOMIC      as GENOMIC      } from './ascc_genomic'
include { ASCC_ORGANELLAR   as ORGANELLAR   } from './ascc_organellar'

include { softwareVersionsToYAML            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText            } from '../subworkflows/local/utils_nfcore_ascc_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow ASCC {

    take:
    genomic_genomes         // channel: samplesheet read in from --input
    organellar_genomes      // channel: tuple(meta, reference)
    fcs_db                  // [path(path)]
    collected_reads         //
    scientific_name         // val(name)
    pacbio_database         // tuple [[meta.id], pacbio_database]
    ncbi_lineage_path

    main:
    ch_versions     = Channel.empty()

    //
    // WORKFLOW: Run main workflow for GENOMIC samples
    //
    GENOMIC (
        genomic_genomes,
        organellar_genomes,
        fcs_db,
        collected_reads,
        scientific_name,
        pacbio_database,
        ncbi_lineage_path
    )
    ch_versions         = ch_versions.mix(GENOMIC.out.versions)


    //
    // WORKFLOW: Run main workflow for ORGANELLAR samples
    //
    if ( params.genomic_only ) {
        ORGANELLAR (
            organellar_genomes,
            fcs_db,
            collected_reads,
            scientific_name,
            pacbio_database,
            ncbi_lineage_path
        )
        ch_versions     = ch_versions.mix(ORGANELLAR.out.versions)
    }


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

}
