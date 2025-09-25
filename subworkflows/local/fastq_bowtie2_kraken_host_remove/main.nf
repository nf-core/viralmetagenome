include { KRAKEN2_KRAKEN2 as KRAKEN2_HOST_REMOVE     } from '../../../modules/nf-core/kraken2/kraken2/main'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_HOST_REMOVE } from '../../../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_HOST_REMOVE } from '../../../modules/nf-core/bowtie2/build/main'
include { FASTQC as FASTQC_HOST                      } from '../../../modules/nf-core/fastqc/main'

workflow FASTQ_BOWTIE2_KRAKEN_HOST_REMOVE {
    take:
    ch_reads
    ch_kraken2_host_db
    host_fna
    tool
    skip_fastqc
    min_reads

    main:
    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()

    if (tool == "kraken2") {
        // remove host reads & keep unclassified reads [true, true]
        KRAKEN2_HOST_REMOVE(
            ch_reads,
            ch_kraken2_host_db,
            true,
            true,
        )

        ch_versions          = ch_versions.mix(KRAKEN2_HOST_REMOVE.out.versions.first())
        ch_multiqc           = ch_multiqc.mix(KRAKEN2_HOST_REMOVE.out.report)
        ch_reads_hostremoved = KRAKEN2_HOST_REMOVE.out.unclassified_reads_fastq
            .join(KRAKEN2_HOST_REMOVE.out.report)
            .map { meta, fastq, tsv -> [meta, fastq, getReadsAfterHostRemoveKraken(tsv)] }
            .branch { meta, fastq, n_reads ->
                pass: n_reads > min_reads
                return [meta, fastq]
                fail: n_reads <= min_reads
                return [meta, n_reads]
            }

    } else {
        // use Bowtie2 to remove
        ch_bowtie2_fasta = host_fna.map{ it -> [ [id:it.baseName], it ]}.collect()
        ch_bowtie2_index = BOWTIE2_BUILD_HOST_REMOVE ( ch_bowtie2_fasta ).index
        ch_versions      = ch_versions.mix( BOWTIE2_BUILD_HOST_REMOVE.out.versions )

        // reads, index, fasta, store_unmapped, sort_bam
        BOWTIE2_ALIGN_HOST_REMOVE(ch_reads, ch_bowtie2_index, ch_bowtie2_fasta, true, false)
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN_HOST_REMOVE.out.versions.first())
        ch_multiqc  = ch_multiqc.mix(BOWTIE2_ALIGN_HOST_REMOVE.out.log)

        ch_reads_hostremoved = BOWTIE2_ALIGN_HOST_REMOVE.out.fastq
            .join(BOWTIE2_ALIGN_HOST_REMOVE.out.log)
            .map { meta, fastq, tsv -> [meta, fastq, getReadsAfterHostRemoveBowtie2(tsv)] }
            .branch { meta, fastq, n_reads ->
                pass: n_reads > min_reads
                return [meta, fastq]
                fail: n_reads <= min_reads
                return [meta, n_reads]
            }
    }

    // fastqc
    if (!skip_fastqc) {
        FASTQC_HOST(
            ch_reads_hostremoved.pass
        )
       ch_multiqc = ch_multiqc.mix(FASTQC_HOST.out.html, FASTQC_HOST.out.zip)
       ch_versions= ch_versions.mix(FASTQC_HOST.out.versions.first())
    }

    emit:
    reads_hostremoved      = ch_reads_hostremoved.pass // channel: [ [ meta ], [ fastq ] ]
    reads_hostremoved_fail = ch_reads_hostremoved.fail // channel: [ [ meta ], [ n_reads ] ]
    mqc                    = ch_multiqc                // channel: [ multiqc_files ]
    versions               = ch_versions               // channel: [ versions.yml ]
}

def getReadsAfterHostRemoveKraken(tsv) {
    def n_reads = 0L
    def firstLine = tsv.readLines().first()
    if (firstLine =~ /(unclassified)/) {
        def valuePart = firstLine.split('\t')[1]
        // Fetching the value from the second column
        n_reads += valuePart.toLong()
    }
    return n_reads.toLong()
}

def getReadsAfterHostRemoveBowtie2(tsv) {
    def n_reads = 0L
    tsv.eachLine { line ->
        def m = line =~ /^\\s*(\\d+) \\([0-9.]+%\\) aligned 0 times/
        if (m.find()) {
            n_reads += m[0][1].toLong()
        }
    }
    return n_reads
}
