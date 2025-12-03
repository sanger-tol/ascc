//
// LOCAL MODULE IMPORT
//
include { DECONTAMINATE_GENERATE_BED        } from '../../../modules/local/decontaminate/generate_bed/main'
include { DECONTAMINATE_CLIP_REGIONS_FASTA  } from '../../../modules/local/decontaminate/clip_regions_fasta/main'

// FUNCTION IMPORTS
// NOTE: IN FUTURE SHOULD ALSO CONTAIN DATA-MAPPER FUNCTIONS
include { getEmptyPlaceholder                           } from "${projectDir}/lib/ascc_utils.groovy"

workflow RUN_DECONTAMINATE_FASTA {
    take:
    input_genome            // Channel.of([ [ id: sample_id ], fasta ])
    parsed_fcsgx
    fcsgx_tiara_summary
    fcs_adaptor_file
    trailingns
    barcodes
    mito_recommendations
    plastid_recommendations

    main:
    ch_versions             = Channel.empty()

    //
    // LOGIC: DATAMAPPER FOR DECONTAMINATE_GENERATE_BED
    //
    input_genome
        .map{meta, file -> [[id: meta.id], file]}
        .join(parsed_fcsgx,             remainder: true)
        .join(fcsgx_tiara_summary,      remainder: true)
        .join(fcs_adaptor_file,         remainder: true)
        .join(trailingns,               remainder: true)
        .join(barcodes,                 remainder: true)
        .join(mito_recommendations,     remainder: true)
        .join(plastid_recommendations,  remainder: true)
        .filter { items ->
            def meta = items[0]
            meta != null &&
            meta != [] &&
            !(meta instanceof Map && (meta.id == null || meta.isEmpty()))
        }
        .map { items ->
            // Replace null values with placeholder file
            items.withIndex().collect { item, index ->
                if (item == null) {
                    getEmptyPlaceholder(index)
                } else if (item instanceof List && item.isEmpty()) {
                    getEmptyPlaceholder(index)
                } else {
                    item
                }
            }
        }
        .set{ merge_input_channel}


    //
    // MODULE: GENERATES CONTAMINATION .BED
    //
    DECONTAMINATE_GENERATE_BED (
        merge_input_channel
    )
    ch_versions             = ch_versions.mix(DECONTAMINATE_GENERATE_BED.out.versions)


    //
    // MODULE: CREATES A FASTA CONTAINING SLIDING WINDOWS OF THE INPUT GENOME
    //
    DECONTAMINATE_CLIP_REGIONS_FASTA (
        DECONTAMINATE_GENERATE_BED.out.main_contamination_data
    )
    ch_versions             = ch_versions.mix(DECONTAMINATE_CLIP_REGIONS_FASTA.out.versions)


    emit:
    // contamination_report_txt    = DECONTAMINATE_CLIP_REGIONS_FASTA.out.contamination_report_txt
    // contamination_bed           = DECONTAMINATE_CLIP_REGIONS_FASTA.out.contamination_bed
    // decontaminated_fasta        = DECONTAMINATE_CLIP_REGIONS_FASTA.out.decontaminated_fasta
    versions                    = ch_versions

}
