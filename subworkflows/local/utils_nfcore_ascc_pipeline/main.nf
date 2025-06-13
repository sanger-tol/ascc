//
// Subworkflow with functionality specific to the sanger-tol/ascc pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

// ASCC Custom param pre-processing
include { VALIDATE_TAXID            } from '../../../modules/local/validate/taxid/main'
include { GUNZIP                    } from '../../../modules/nf-core/gunzip/main'
include { PREPARE_BLASTDB           } from '../../local/prepare_blastdb/main'
include { CHECK_NT_BLAST_TAXONOMY   } from '../../../modules/local/check/nt_blast_taxonomy/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input
    //

    Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map {
            sample, type_of_assembly, assembly_file ->
                return [[ id: sample.id + '_' + type_of_assembly, assembly_type: type_of_assembly, single_end:true ], assembly_file ]
        }
        .groupTuple()
        // .map { samplesheet ->
        //     validateInputSamplesheet(samplesheet)
        // }
        .map {
            meta, fastas ->
                return [ meta, fastas[0] ] // We are only expecting one fasta file per sample+haplo
        }
        .set { ch_samplesheet }


    //
    // MODULE: ENSURE THAT THE TAXID FOR THE INPUT GENOME IS INDEED IN THE TAXDUMP
    //         MODULE DOES NOT OUTPUT ANYTHING BUT SHOULD KILL PIPELINE ON FAIL
    //
    VALIDATE_TAXID(
        Channel.of(params.taxid),
        Channel.of(params.ncbi_taxonomy_path)
    )
    ch_versions = ch_versions.mix(VALIDATE_TAXID.out.versions)


    //
    // LOGIC: GUNZIP INPUT DATA IF GZIPPED, OTHERWISE PASS
    //
    ch_samplesheet
        .branch { meta, file ->
            zipped: file.name.endsWith('.gz')
            unzipped: !file.name.endsWith('.gz')
        }
        .set {ch_input}


    //
    // MODULE: UNZIP INPUTS IF NEEDED
    //
    GUNZIP (
        ch_input.zipped
    )
    ch_versions = ch_versions.mix(GUNZIP.out.versions)


    //
    // LOGIC: MIX CHANELS WHICH MAY OR MAY NOT BE EMPTY INTO A SINGLE QUEUE CHANNEL
    //
    unzipped_input = Channel.empty()

    unzipped_input
        .mix(ch_input.unzipped, GUNZIP.out.gunzip)
        .set { standardised_unzipped_input }


    //
    // LOGIC: FILTER THE INPUT BASED ON THE assembly_type VALUE IN THE META
    //          DEPENDING ON THIS VALUE THE PIPELINE WILL NEED TO BE DIFFERENT
    //
    standardised_unzipped_input
        .branch{
            organellar_genome: it[0].assembly_type == "MITO" || it[0].assembly_type == "PLASTID"
            genomic_genome: it[0].assembly_type  == "PRIMARY" || it[0].assembly_type  == "HAPLO"
            error: true
        }
        .set { branched_assemblies }


    //
    // NOTE: Setting the basic channels form the input
    //
    Channel.fromPath(params.pacbio_barcode_file)
        .set {barcode_data_file}

    Channel.of(params.fcs_gx_database_path)
        .set { fcs_gx_database_path}




    //
    // SUBWORKFLOW: PREPARE THE MAKEBLASTDB INPUTS
    //
    PREPARE_BLASTDB (
        params.sample_id,
        params.reads_path,
        params.reads_type,
        barcode_data_file,
        params.pacbio_barcode_names
    )
    versions = ch_versions.mix(PREPARE_BLASTDB.out.versions)


    //
    // LOGIC: GETS PACBIO READ PATHS FROM READS_PATH IF (COVERAGE OR BTK SUBWORKFLOW IS ACTIVE) OR ALL
    //
    if ( params.run_coverage != "off" || params.run_btk != "off" ) {
        ch_grabbed_reads_path       = Channel.of(params.reads_path).collect()
    } else {
        ch_grabbed_reads_path       = []
    }


    //
    // MODULE: CHECK IF NT BLAST DATABASE HAS TAXONOMY INCLUDED (ONLY IF NT BLAST IS INCLUDED)
    // This check is specifically for the nt BLAST database used in the EXTRACT_NT_BLAST subworkflow,
    // not for other BLAST databases used elsewhere in the pipeline (VecScreen, PacBio barcodes check, etc.)
    //
    if ( params.run_nt_blast != "off" ) {
        CHECK_NT_BLAST_TAXONOMY(
            params.nt_database_path
        )
        ch_versions = ch_versions.mix(CHECK_NT_BLAST_TAXONOMY.out.versions)

        // Check the result and fail if needed
        CHECK_NT_BLAST_TAXONOMY.out.status
            .map { it.trim() }  // Trim any whitespace
            .subscribe { status ->
                if (status == "nt_database_taxonomy_files_not_found") {
                    log.error "NT BLAST database taxonomy check failed"
                    exit 1, "The NT BLAST database does not have taxonomy included. Please see the error message above for details."
                }
            }
    }


    emit:
    samplesheet             = ch_samplesheet
    sample_id               = params.sample_id
    main_genomes            = branched_assemblies.genomic_genome
    organellar_genomes      = branched_assemblies.organellar_genome
    barcodes_file           = barcode_data_file
    pacbio_db               = PREPARE_BLASTDB.out.barcodes_blast_db
    fcs_gx_database         = fcs_gx_database_path
    collected_reads         = ch_grabbed_reads_path
    versions                = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                []
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate channels from input samplesheet
//
// def validateInputSamplesheet(input) {
//     def (metas, fastqs) = input[1..2]

//     // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
//     def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
//     if (!endedness_ok) {
//         error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
//     }

//     return [ metas[0], fastqs ]
// }

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
