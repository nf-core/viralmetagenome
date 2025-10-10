//
// Consensus calling with BCFTools
//

include { TABIX_TABIX        } from '../../../modules/nf-core/tabix/tabix/main'
include { BEDTOOLS_MERGE     } from '../../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_MASKFASTA } from '../../../modules/nf-core/bedtools/maskfasta/main'
include { BCFTOOLS_CONSENSUS } from '../../../modules/nf-core/bcftools/consensus/main'
include { MAKE_BED_MASK      } from '../../../modules/local/make_bed_mask/main'

workflow BAM_VCF_CONSENSUS_BCFTOOLS {
    take:
    ch_bam        // channel: [ val(meta), [ bam ] ]
    ch_vcf        // channel: [ val(meta), [ vcf ] ]
    ch_fasta      // channel: [ val(meta), [ fasta ] ]
    mapping_stats // value: [ true | false ]

    main:

    ch_versions = Channel.empty()

    TABIX_TABIX(
        ch_vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    ch_bam_vcf_fasta = ch_bam
        .join(ch_vcf, by: [0])
        .join(ch_fasta, by: [0])

    //
    // Create BED file with consensus regions to mask (regions to remove)
    //
    MAKE_BED_MASK(
        ch_bam_vcf_fasta,
        mapping_stats,
    )
    ch_versions = ch_versions.mix(MAKE_BED_MASK.out.versions.first())

    //
    // Merge intervals with BEDTools
    //
    BEDTOOLS_MERGE(
        MAKE_BED_MASK.out.bed
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE.out.versions.first())

    ch_bed_fasta = BEDTOOLS_MERGE.out.bed.join(ch_fasta, by: [0])

    //
    // Mask regions in consensus with BEDTools
    //
    BEDTOOLS_MASKFASTA(
        ch_bed_fasta
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MASKFASTA.out.versions.first())

    //
    // Call consensus sequence with BCFTools
    //
    ch_bcftools_in = ch_vcf
        .join(TABIX_TABIX.out.tbi, by: [0])
        .join(BEDTOOLS_MASKFASTA.out.fasta, by: [0])
        .map { meta, v, t, f ->
            [meta, v, t, f, []]
        }
    BCFTOOLS_CONSENSUS(
        ch_bcftools_in
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions.first())

    emit:
    consensus = BCFTOOLS_CONSENSUS.out.fasta // channel: [ val(meta), [ fasta ] ]
    versions  = ch_versions // channel: [ versions.yml ]
}
