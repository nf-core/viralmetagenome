// Based on https://github.com/nf-core/viralrecon/blob/master/subworkflows/local/variants_ivar.nf

include { SNPEFF_BUILD          } from '../../modules/local/snpeff_build/main'
include { SNPEFF_SNPEFF         } from '../../modules/nf-core/snpeff/snpeff/main'
include { SNPSIFT_EXTRACTFIELDS } from '../../modules/local/snpsift_extractfields/main'

workflow VCF_ANNOTATE {

    take:
    ch_vcf_ref         // channel: [ val(meta), [ vcf], [ref]]

    main:
    ch_versions = Channel.empty()

    // A gff must be supplied & the genome shouldn't have undergone selection
    ch_vcf_ref = ch_vcf_ref
        .filter{ meta, vcf, ref -> meta.gff != null && !meta.selection  }

    ch_gff = ch_vcf_ref
        .map{ meta, vcf, ref -> [ ref, meta.gff ] }
        .unique()

    SNPEFF_BUILD (
        ch_gff
    )
    ch_versions = ch_versions.mix(SNPEFF_BUILD.out.versions.first())

    ch_vcf_ref
        .map{ meta, vcf, ref -> [ ref, meta, vcf ] }
        .combine(SNPEFF_BUILD.out.db, by: [0])
        .multiMap { ref, meta, vcf, db_id, db, config ->
            vcf: [meta, vcf]
            db: db_id
            cache: [meta, db, config]
        }
        .set { snpeff_in }

    SNPEFF_SNPEFF (
        snpeff_in.vcf,
        snpeff_in.db,
        snpeff_in.cache,
    )
    ch_vcf_ann  = SNPEFF_SNPEFF.out.vcf
    ch_versions = ch_versions.mix(SNPEFF_SNPEFF.out.versions.first())

    SNPSIFT_EXTRACTFIELDS (
        ch_vcf_ann,
    )
    ch_versions = ch_versions.mix(SNPSIFT_EXTRACTFIELDS.out.versions.first())


    emit:
    vcf          = ch_vcf_ann                     // channel: [ val(meta), [ tsv ] ]
    txt          = SNPSIFT_EXTRACTFIELDS.out.txt  // channel: [ val(meta), [ txt ] ]

    versions     = ch_versions                    // channel: [ versions.yml ]
}

