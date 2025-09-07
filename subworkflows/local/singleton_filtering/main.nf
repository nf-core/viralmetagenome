include { filterContigs                                        } from '../utils_nfcore_viralmetagenome_pipeline'
include { RENAME_FASTA_HEADER as RENAME_FASTA_HEADER_SINGLETON } from '../../../modules/local/rename_fasta_header/main'

workflow SINGLETON_FILTERING {

    take:
    ch_fasta               // channel: [ val(meta), [ fasta ] ]
    min_contig_size     // int
    max_n_perc        // int

    main:
    ch_versions = Channel.empty()

    if ( !params.skip_singleton_filtering) {
        filtered = filterContigs ( ch_fasta, min_contig_size, max_n_perc)
        ch_contig = filtered.pass
    }
    // Rename to avoid errors downstream
    RENAME_FASTA_HEADER_SINGLETON(
        ch_contig,
        []
        )
    ch_versions = ch_versions.mix(RENAME_FASTA_HEADER_SINGLETON.out.versions)


    emit:
    filtered     = RENAME_FASTA_HEADER_SINGLETON.out.fasta  // channel: [ val(meta), [ fasta ] ]
    versions     = ch_versions                              // channel: [ versions.yml ]
}
