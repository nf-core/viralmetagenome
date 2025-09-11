include { UNTAR as UNTAR_DB      } from '../../../modules/nf-core/untar/main'
include { GUNZIP as GUNZIP_DB    } from '../../../modules/nf-core/gunzip/main'
include { XZ_DECOMPRESS as XZ_DB } from '../../../modules/nf-core/xz/decompress/main'

workflow UNPACK_DB {
    take:
    ch_db_compressed // channel [ val(meta), [ db ] ]

    main:
    ch_versions = Channel.empty()
    ch_db_unpacked = Channel.empty()

    ch_db_branched = ch_db_compressed.branch { _meta, dbs ->
        tar: dbs.name.endsWith('.tar.gz') || dbs.name.endsWith('.tgz') || dbs.name.endsWith('.tar')
        gzip: dbs.name.endsWith('.gz')
        xz: dbs.name.endsWith('.xz')
        other: true
    }

    ch_db_unpacked = ch_db_unpacked.mix(ch_db_branched.other)

    ch_db_unpacked = ch_db_unpacked.mix(UNTAR_DB(ch_db_branched.tar).untar)
    ch_versions = ch_versions.mix(UNTAR_DB.out.versions.first())

    ch_db_unpacked = ch_db_unpacked.mix(GUNZIP_DB(ch_db_branched.gzip).gunzip)
    ch_versions = ch_versions.mix(GUNZIP_DB.out.versions.first())

    ch_db_unpacked = ch_db_unpacked.mix(XZ_DB(ch_db_branched.xz).file)
    ch_versions = ch_versions.mix(XZ_DB.out.versions.first())

    emit:
    db       = ch_db_unpacked // channel: [ db ]
    versions = ch_versions    // channel: [ versions.yml ]
}
