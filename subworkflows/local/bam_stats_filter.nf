// Take in a bam file and remove those that don't have any reads aligned

include {failedMappedReadsToMultiQC } from '../../subworkflows/local/utils_nfcore_viralgenie_pipeline'
include { SAMTOOLS_INDEX            } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS            } from '../../modules/nf-core/samtools/stats/main'

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

    // Create clean copies before joining
    ch_bam
        .map { meta, bam -> [meta + [:], bam] }
        .set { ch_bam_clean }

    SAMTOOLS_INDEX.out.bai
        .map { meta, bai -> [meta + [:], bai] }
        .set { ch_bai_clean }

    ch_reference
        .map { meta, ref -> [meta + [:], ref] }
        .set { ch_reference_clean }

    stats_in = ch_bam_clean
        .join(ch_bai_clean, by: [0])
        .join(ch_reference_clean, by: [0])
        .multiMap{ meta, bam, bai, ref ->
            bam_bai : [meta, bam, bai ]
            ref : [meta, ref ]
        }

    SAMTOOLS_STATS ( stats_in.bam_bai, stats_in.ref )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())

    // Create clean copies before joining
    SAMTOOLS_STATS
        .out
        .stats
        .map { meta, stats -> [meta + [:], stats] }
        .set { ch_stats_clean }

    ch_bam_clean
        .join(ch_stats_clean, by: [0] )
        .map{ meta, bam, stats -> [ meta, bam, WorkflowCommons.getStatsMappedReads(stats) ] }
        .branch { meta, bam, mapped_reads ->
            pass: mapped_reads > min_mapped_reads
                return [ meta, bam ]
            fail: mapped_reads <= min_mapped_reads
                return [ meta, bam, mapped_reads ]
        }
        .set{ ch_bam_filtered }

    bam_pass = ch_bam_filtered.pass
    bam_fail = ch_bam_filtered.fail

    ch_fail_mapping_multiqc = failedMappedReadsToMultiQC(bam_fail, min_mapped_reads).collectFile(name:'failed_mapped_reads_mqc.tsv')

    emit:
    bam_pass     = bam_pass                    // channel: [ val(meta), [ bam ] ]
    stats        = SAMTOOLS_STATS.out.stats    // channel: [ val(meta), [ stats ] ]
    bam_fail_mqc = ch_fail_mapping_multiqc     // channel: [ val(meta), [ bam ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

