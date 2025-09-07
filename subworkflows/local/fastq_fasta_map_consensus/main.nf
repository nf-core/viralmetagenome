include { filterContigs; failedContigsToMultiQC   } from '../utils_nfcore_viralmetagenome_pipeline'
include { MAP_READS                               } from '../map_reads'
include { BAM_DEDUPLICATE                         } from '../bam_deduplicate'
include { SAMTOOLS_FAIDX                          } from '../../../modules/nf-core/samtools/faidx/main'
include { BAM_STATS_METRICS                       } from '../bam_stats_metrics'
include { BAM_CALL_VARIANTS                       } from '../bam_call_variants'
include { BAM_CALL_CONSENSUS                      } from '../bam_call_consensus'
include { BAM_STATS_FILTER                        } from '../bam_stats_filter'

workflow FASTQ_FASTA_MAP_CONSENSUS {

    take:
    ch_reference_reads   // channel: [ val(meta), [ fasta ], [ fastq ] ]
    mapper               // val: [ bwamem2 | bowtie2 ]
    umi                  // val: [ true | false ]
    deduplicate          // val: [ true | false ]
    call_variants        // val: [ true | false ]
    variant_caller       // val: [ bcftools | ivar ]
    consensus_caller     // val: [ bcftools | ivar ]
    mapping_stats        // val: [ true | false ]
    min_mapped_reads     // integer: min_mapped_reads
    min_len              // integer: min_length
    n_100                // integer: n_100

    main:

    ch_versions     = Channel.empty()
    ch_multiqc      = Channel.empty()
    ch_dedup_bam    = Channel.empty()
    ch_reads_in     = ch_reference_reads.map{meta, ref, reads -> [meta,reads] }

    // mapping of reads using bowtie2 or BWA-MEM2
    MAP_READS ( ch_reference_reads, mapper )

    ch_bam       = MAP_READS.out.bam
    ch_reference = MAP_READS.out.ref
    ch_versions  = ch_versions.mix(MAP_READS.out.versions)
    ch_multiqc   = ch_multiqc.mix(MAP_READS.out.mqc.collect{it[1]}.ifEmpty([]))

    SAMTOOLS_FAIDX ( ch_reference, [[],[]], false)
    ch_versions  = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    // remove references-read combinations with low mapping rates
    BAM_STATS_FILTER ( ch_bam, ch_reference, min_mapped_reads )
    ch_multiqc   = ch_multiqc.mix(BAM_STATS_FILTER.out.stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc   = ch_multiqc.mix(BAM_STATS_FILTER.out.bam_fail_mqc.ifEmpty([]))
    ch_versions  = ch_versions.mix(BAM_STATS_FILTER.out.versions)

     ch_bam_fa_fai = BAM_STATS_FILTER.out.bam_pass
        .join(ch_reference, by: [0])
        .join(SAMTOOLS_FAIDX.out.fai, by: [0])

    // deduplicate bam using umitools (if UMI) or picard
    if (deduplicate) {
        BAM_DEDUPLICATE ( ch_bam_fa_fai, umi, mapping_stats)

        ch_dedup_bam = BAM_DEDUPLICATE.out.bam
        ch_multiqc   = ch_multiqc.mix(BAM_DEDUPLICATE.out.mqc.collect{it[1]}.ifEmpty([]))
        ch_versions  = ch_versions.mix(BAM_DEDUPLICATE.out.versions)

    } else {
        ch_dedup_bam = ch_bam
    }

    ch_dedup_bam_ref = ch_dedup_bam
        .join(ch_reference, by: [0])

    // report summary statistics of alignment
    if (mapping_stats) {
        BAM_STATS_METRICS ( ch_dedup_bam_ref )
        ch_multiqc   = ch_multiqc.mix(BAM_STATS_METRICS.out.mqc.collect{it[1]}.ifEmpty([]))
        ch_versions  = ch_versions.mix(BAM_STATS_METRICS.out.versions)
    }

    // call variants
    ch_vcf        = Channel.empty()
    ch_vcf_filter = Channel.empty()
    ch_tbi        = Channel.empty()

    if (consensus_caller == "bcftools" || call_variants ) {
        BAM_CALL_VARIANTS (
            ch_dedup_bam_ref,
            variant_caller,
            mapping_stats
        )
        ch_versions   = ch_versions.mix(BAM_CALL_VARIANTS.out.versions)
        ch_multiqc    = ch_multiqc.mix(BAM_CALL_VARIANTS.out.mqc.collect{it[1]}.ifEmpty([]))
        ch_vcf_filter = BAM_CALL_VARIANTS.out.vcf_filter
        ch_vcf        = BAM_CALL_VARIANTS.out.vcf
        ch_tbi        = BAM_CALL_VARIANTS.out.tbi
    }

    // cannot merge ch_vcf_filter as it will not have a meta

    // consensus calling
    BAM_CALL_CONSENSUS (
        ch_dedup_bam_ref,
        ch_vcf_filter,
        consensus_caller,
        mapping_stats
    )
    ch_versions = ch_versions.mix(BAM_CALL_CONSENSUS.out.versions)
    ch_consensus_all      = BAM_CALL_CONSENSUS.out.consensus

    // Check if consensus genomes are long enough
    ch_contigs            = filterContigs ( ch_consensus_all, min_len, n_100 )
    ch_consensus_filtered = ch_contigs.pass
    ch_contig_qc_fail_mqc = failedContigsToMultiQC ( ch_contigs.fail, min_len, n_100 )
    ch_consensus_reads    = ch_consensus_filtered.join(ch_reads_in, by: [0])

    // Vcf & bam files
    ch_vcf_ref         = ch_vcf.join(ch_reference, by: [0])
    bam_out            = ch_dedup_bam_ref.map{meta,bam,ref -> [meta,bam] }

    ch_multiqc         = ch_multiqc.mix(ch_contig_qc_fail_mqc.collectFile(name:'failed_contig_quality_mqc.tsv').ifEmpty([]))

    emit:
    ch_consensus_reads = ch_consensus_reads                     // channel: [ val(meta), [ fasta ], [ fastq ] ]
    ch_consensus       = ch_consensus_filtered                  // channel: [ val(meta), [ fasta ] ]
    ch_consensus_all   = ch_consensus_all                       // channel: [ val(meta), [ fasta ] ]
    ch_bam             = bam_out                                // channel: [ val(meta), [ bam ] ]
    ch_vcf             = ch_vcf                                 // channel: [ val(meta), [ vcf ] ]
    ch_vcf_ref         = ch_vcf_ref                             // channel: [ val(meta), [ vcf ], [ fasta ] ]
    ch_vcf_filter      = ch_vcf_filter                          // channel: [ val(meta), [ vcf ] ]

    ch_mqc             = ch_multiqc                             // channel: [ val(meta), [ csi ] ]
    ch_versions        = ch_versions                            // channel: [ versions.yml ]
}
