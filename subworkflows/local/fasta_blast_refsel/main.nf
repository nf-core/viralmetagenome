include { BLAST_BLASTN          } from '../../../modules/nf-core/blast/blastn/main'
include { BLAST_FILTER          } from '../../../modules/local/blast_filter'

workflow FASTA_BLAST_REFSEL {

    take:
    ch_fasta          // channel: [ val(meta), path(fasta)]
    ch_blast_db       // channel: [ val(meta), path(db) ]
    ch_blast_db_fasta // channel: [ val(meta), path(fasta) ]

    main:
    ch_versions = Channel.empty()
    // Blast results, to a reference database, to find a complete genome that's already assembled
    BLAST_BLASTN (
        ch_fasta,
        ch_blast_db
    )
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())

    ch_blast_txt = BLAST_BLASTN.out.txt
        .branch { _meta, txt ->
            no_hits : txt.countLines() == 0
            hits    : txt.countLines() > 0
        }

    // Make a table of samples that did not have any blast hits
    ch_no_blast_hits = Channel.empty()
    ch_no_blast_hits = ch_blast_txt
        .no_hits
        .join(ch_fasta)


    // Filter out false positve hits that based on query length, alignment length, identity, e-score & bit-score
    ch_blast_txt
        .hits
        .join(ch_fasta, by:[0], remainder:true)
        .multiMap{
            meta, txt, fasta ->
            hits : [meta, txt ? txt : []]
            contigs : [meta, fasta]
        }
        .set{input_blast_filter}

    BLAST_FILTER (
        input_blast_filter.hits,
        input_blast_filter.contigs,
        ch_blast_db_fasta
    )
    ch_versions = ch_versions.mix(BLAST_FILTER.out.versions.first())

    emit:
    fasta_ref_contigs = BLAST_FILTER.out.sequence  // channel: [ val(meta), [ fasta ] ]
    no_blast_hits     = ch_no_blast_hits           // channel: [ val(meta), [ mqc ] ]
    versions          = ch_versions                // channel: [ versions.yml ]
}
