//
// Consensus calling with BCFTools
//

include { TABIX_TABIX         } from '../../modules/nf-core/tabix/tabix/main'
include { BEDTOOLS_MERGE      } from '../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_MASKFASTA  } from '../../modules/nf-core/bedtools/maskfasta/main'
include { BCFTOOLS_CONSENSUS  } from '../../modules/nf-core/bcftools/consensus/main'
include { MAKE_BED_MASK       } from '../../modules/local/make_bed_mask/main'

workflow BAM_VCF_CONSENSUS_BCFTOOLS {
    take:
    ch_bam          // channel: [ val(meta), [ bam ] ]
    ch_vcf          // channel: [ val(meta), [ vcf ] ]
    ch_fasta        // channel: [ val(meta), [ fasta ] ]
    mapping_stats    // value: [ true | false ]

    main:

    ch_versions = Channel.empty()

    TABIX_TABIX (
        ch_vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    // Create clean copies before joining
    ch_bam
        .map { meta, bam -> [meta.clone(), bam] }
        .set { bam_clean }

    ch_vcf
        .map { meta, vcf -> [meta.clone(), vcf] }
        .set { vcf_clean }

    ch_fasta
        .map { meta, fasta -> [meta.clone(), fasta] }
        .set { fasta_clean }

    bam_clean
        .join(vcf_clean, by: [0])
        .join(fasta_clean, by: [0])
        .set{bam_vcf_fasta}

    //
    // Create BED file with consensus regions to mask (regions to remove)
    //
    MAKE_BED_MASK (
        bam_vcf_fasta,
        mapping_stats
    )
    ch_versions = ch_versions.mix(MAKE_BED_MASK.out.versions.first())

    //
    // Merge intervals with BEDTools
    //
    BEDTOOLS_MERGE (
        MAKE_BED_MASK.out.bed
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE.out.versions.first())

    // Create clean copies before joining
    BEDTOOLS_MERGE
        .out
        .bed
        .map { meta, bed -> [meta.clone(), bed] }
        .join(fasta_clean, by: [0])
        .set{bed_fasta}

    //
    // Mask regions in consensus with BEDTools
    //
    BEDTOOLS_MASKFASTA (
        bed_fasta
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MASKFASTA.out.versions.first())

    //
    // Call consensus sequence with BCFTools
    //
    // Create clean copies before joining
    vcf_clean
        .map { meta, vcf -> [meta.clone(), vcf] }
        .set { vcf_join_clean }

    TABIX_TABIX.out.tbi
        .map { meta, tbi -> [meta.clone(), tbi] }
        .set { tbi_clean }

    BEDTOOLS_MASKFASTA.out.fasta
        .map { meta, fasta -> [meta.clone(), fasta] }
        .set { masked_fasta_clean }

    BCFTOOLS_CONSENSUS (
        vcf_join_clean.join(tbi_clean, by: [0]).join(masked_fasta_clean, by: [0])
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions.first())

    emit:
    consensus        = BCFTOOLS_CONSENSUS.out.fasta     // channel: [ val(meta), [ fasta ] ]
    versions         = ch_versions                       // channel: [ versions.yml ]
}
