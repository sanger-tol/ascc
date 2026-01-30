/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ASCC_GENOMIC      as GENOMIC      } from './ascc_genomic'
include { ASCC_ORGANELLAR   as ORGANELLAR   } from './ascc_organellar'

include { softwareVersionsToYAML            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText            } from '../subworkflows/local/utils_nfcore_ascc_pipeline'
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
    val_reads_per_chunk

    main:
    ch_versions     = channel.empty()

    //
    // WORKFLOW: Run main workflow for GENOMIC samples
    //
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
        barcodes,
        val_reads_per_chunk
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
            barcodes,
            val_reads_per_chunk
        )
        ch_versions     = ch_versions.mix(ORGANELLAR.out.versions)
    }


    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'ascc_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_collated_versions                 // channel: [ path(versions.yml) ]

}
