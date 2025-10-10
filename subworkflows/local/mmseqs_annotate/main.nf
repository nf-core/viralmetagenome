include { MMSEQS_CREATEDB as MMSEQS_CREATEANNOTATIONDB } from '../../../modules/nf-core/mmseqs/createdb/main'
include { MMSEQS_EASYSEARCH                            } from '../../../modules/nf-core/mmseqs/easysearch/main'

workflow MMSEQS_ANNOTATE {
    take:
    ch_genomes // channel: [ val(meta), [ fasta ] ]
    ch_db      // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    // create mmseqs annotation db
    MMSEQS_CREATEANNOTATIONDB(ch_db)
    ch_versions = ch_versions.mix(MMSEQS_CREATEANNOTATIONDB.out.versions.first())

    // search the genomes against the annotation db
    MMSEQS_EASYSEARCH(ch_genomes, MMSEQS_CREATEANNOTATIONDB.out.db)
    ch_versions = ch_versions.mix(MMSEQS_EASYSEARCH.out.versions.first())

    emit:
    tsv      = MMSEQS_EASYSEARCH.out.tsv // channel: [ val(meta), [ tsv ] ]
    versions = ch_versions               // channel: [ versions.yml ]
}
