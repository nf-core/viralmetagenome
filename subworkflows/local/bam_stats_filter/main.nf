// Take in a bam file and remove those that don't have any reads aligned

include { failedMappedReadsToMultiQC } from '../utils_nfcore_viralmetagenome_pipeline'
include { getStatsMappedReads        } from '../utils_nfcore_viralmetagenome_pipeline'
include { SAMTOOLS_INDEX             } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS             } from '../../../modules/nf-core/samtools/stats/main'

workflow BAM_STATS_FILTER {

    take:
    ch_bam           // channel: [ val(meta), [ bam ] ]
    ch_reference     // channel: [ val(meta), [ fasta ] ]
    min_mapped_reads // integer: min_mapped_reads

    main:

    ch_versions             = Channel.empty()
    ch_fail_mapping_multiqc = Channel.empty()

    SAMTOOLS_INDEX ( ch_bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    ch_stats_in = ch_bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0])
        .join(ch_reference, by: [0])
        .multiMap{ meta, bam, bai, ref ->
            bam_bai : [meta, bam, bai ]
            ref : [meta, ref ]
        }

    SAMTOOLS_STATS ( ch_stats_in.bam_bai, ch_stats_in.ref )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())

    ch_bam_filtered = SAMTOOLS_STATS
        .out
        .stats
        .join(ch_bam, by: [0] )
        .map{ meta, stats, bam -> [ meta, bam, getStatsMappedReads(stats) ] }
        .branch { meta, bam, mapped_reads ->
            pass: mapped_reads > min_mapped_reads
                return [ meta, bam ]
            fail: mapped_reads <= min_mapped_reads
                return [ meta, bam, mapped_reads ]
        }

    bam_pass = ch_bam_filtered.pass
    bam_fail = ch_bam_filtered.fail

    ch_fail_mapping_multiqc = failedMappedReadsToMultiQC(bam_fail, min_mapped_reads).collectFile(name:'failed_mapped_reads_mqc.tsv')

    emit:
    bam_pass     = bam_pass                    // channel: [ val(meta), [ bam ] ]
    stats        = SAMTOOLS_STATS.out.stats    // channel: [ val(meta), [ stats ] ]
    bam_fail_mqc = ch_fail_mapping_multiqc     // channel: [ val(meta), [ bam ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}
