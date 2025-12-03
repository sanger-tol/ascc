//
// LOCAL MODULE IMPORT
//
include { DECONTAMINATE_CLIP_REGIONS_FASTA  } from '../../../modules/local/decontaminate/clip_regions_fasta/main'
include { DECONTAMINATE_GENERATE_BED        } from '../../../modules/local/decontaminate/generate_bed/main'


workflow DECONTAMINATE_FASTA {
    take:
    input_genome            // Channel.of([ [ id: sample_id ], fasta ])


    main:
    ch_versions             = Channel.empty()

    //
    // LOGIC: DATAMAPPER FOR DECONTAMINATE_GENERATE_BED
    //


    //
    // MODULE: GENERATES CONTAMINATION .BED
    //
    DECONTAMINATE_GENERATE_BED (
        input_genome
    )
    ch_versions             = ch_versions.mix(SEQKIT_SLIDING.out.versions)


    //
    // MODULE: CREATES A FASTA CONTAINING SLIDING WINDOWS OF THE INPUT GENOME
    //
    DECONTAMINATE_CLIP_REGIONS_FASTA ( input_genome )
    ch_versions             = ch_versions.mix(SEQKIT_SLIDING.out.versions)


    emit:
    contamination_report_txt    = DECONTAMINATE_CLIP_REGIONS_FASTA.out.contamination_report_txt
    contamination_bed           = DECONTAMINATE_CLIP_REGIONS_FASTA.out.contamination_bed
    decontaminated_fasta        = DECONTAMINATE_CLIP_REGIONS_FASTA.out.decontaminated_fasta
    versions                    = ch_versions

}
