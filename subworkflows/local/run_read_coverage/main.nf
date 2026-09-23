//
// LOCAL SUBWORKFLOW IMPORTS
//
include { FASTX_MAP_LONG_READS as SE_MAPPING            } from '../../sanger-tol/fastx_map_long_reads/main'
include { PE_MAPPING                                    } from '../pe_mapping/main'

//
// NF-CORE MODULE IMPORTS
//
include { SAMTOOLS_SORT                                 } from '../../../modules/nf-core/samtools/sort/main'
include { COVERM_CONTIG                                 } from '../../../modules/nf-core/coverm/contig/main'

workflow RUN_READ_COVERAGE {

    take:
    reference_tuple          // Channel [ val(meta), path(file) ]
    reads
    platform                 // Channel val( str )
    val_reads_per_chunk

    main:
    ch_versions     = channel.empty()
    ch_align_bam    = channel.empty()
    ch_out_bam      = channel.empty()

    //
    // LOGIC: GETS PACBIO READ PATHS FROM READS_PATH
    //
    collection_of_reads = reads.flatten()

    ref_and_data        =   reference_tuple
                                .combine(collection_of_reads)

    //
    // LOGIC: CHECK IF THE INPUT READ FILE IS PAIRED END OR SINGLE END BASED ON THE READ PLATFORM
    // THEN RUN MINIMAP
    // - Removed the mix function from this as it is not needed, there shouldn't be multiple read
    // types
    //
    //

    if ( platform in ["hifi", "clr", "ont"] ) {

        //
        // MODULE: RUN SINGLE END MAPPING ON THE REFERENCE AND LONGREAD DATA
        //
        ch_map_long_reads_input = ref_and_data
            .groupTuple(by: [0, 1]) // the reads are not a list so we get multiple input channels otherwise
            .multiMap { meta, reference, read_files ->
                def meta_new = meta + [readtype: platform]
                reference: [ meta_new, reference ]
                read_ch: [ meta_new, read_files ]
            }

        SE_MAPPING(
            ch_map_long_reads_input.reference,
            ch_map_long_reads_input.read_ch,
            val_reads_per_chunk,
            true
        )
        ch_versions = ch_versions.mix(SE_MAPPING.out.versions)
        ch_out_bam  = ch_out_bam.mix(SE_MAPPING.out.bam)

    } else if ( platform in ["illumina"] ) {

        //
        // MODULE: RUN PAIRED END MAPPING ON THE REFERENCE AND SHORTREAD DATA
        //
        PE_MAPPING  (
            ref_and_data
        )
        ch_versions     = ch_versions.mix(PE_MAPPING.out.versions)
        ch_align_bam    = PE_MAPPING.out.mapped_bam

        //
        // MODULE: SORT THE MAPPED BAM
        //
        SAMTOOLS_SORT (
            ch_align_bam,
            [[],[]],
            "csi"
        )
        ch_versions = ch_versions.mix( SAMTOOLS_SORT.out.versions )
        ch_out_bam  = ch_out_bam.mix(SAMTOOLS_SORT.out.bam)
    }

    //
    // MODULE: GET READ DEPTH ACROSS THE GENOME
    //
    COVERM_CONTIG (
        ch_out_bam,
        [[],[]],
        true,
        false
    )
    ch_versions = ch_versions.mix(COVERM_CONTIG.out.versions)

    tsv_ch              = COVERM_CONTIG.out.coverage
                            .map { meta, file -> [ [id: meta.id] , file ] }
                            .ifEmpty { [[:],[]] }

    bam_ch              = ch_out_bam
                            .map { meta, file -> [ [id: meta.id] , file ] }
                            .ifEmpty { [[:],[]] }

    emit:
    tsv_ch
    bam_ch
    versions            = ch_versions
}
