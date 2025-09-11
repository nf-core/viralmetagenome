include { BWAMEM2_MEM   } from '../../../modules/nf-core/bwamem2/mem/main'
include { BWAMEM2_INDEX } from '../../../modules/nf-core/bwamem2/index/main'
include { BOWTIE2_ALIGN } from '../../../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_BUILD } from '../../../modules/nf-core/bowtie2/build/main'
include { BWA_MEM       } from '../../../modules/nf-core/bwa/mem/main'
include { BWA_INDEX     } from '../../../modules/nf-core/bwa/index/main'

workflow MAP_READS {
    take:
    ch_reference_reads // channel: [ val(meta), [ fasta ], [ reads ] ]
    mapper             // val: 'bwamem2' or 'bowtie2' or 'bwa'

    main:

    ch_versions = Channel.empty()
    ch_multiqc = Channel.empty()

    ch_reads = ch_reference_reads.map { meta, fasta, fastq -> [meta, fastq] }
    ch_reference = ch_reference_reads.map { meta, fasta, fastq -> [meta, fasta] }

    if (mapper == 'bwamem2') {
        BWAMEM2_INDEX(ch_reference)
        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions.first())

        ch_bwamem2_input = ch_reference_reads
            .join(BWAMEM2_INDEX.out.index, by: [0])
            .multiMap { meta, fasta, fastq, index ->
                reads: [meta, fastq]
                index: [meta, index]
                fasta: [meta, fasta]
            }

        BWAMEM2_MEM(ch_bwamem2_input.reads, ch_bwamem2_input.index, ch_bwamem2_input.fasta, true)
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())
        //no mqc for bwamem2

        ch_bam = BWAMEM2_MEM.out.bam
    }
    else if (mapper == 'bowtie2') {
        BOWTIE2_BUILD(ch_reference)
        ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions.first())

        ch_bowtie2_input = ch_reference_reads
            .join(BOWTIE2_BUILD.out.index, by: [0])
            .multiMap { meta, fasta, fastq, index ->
                reads: [meta, fastq]
                index: [meta, index]
                fasta: [meta, fasta]
            }

        BOWTIE2_ALIGN(ch_bowtie2_input.reads, ch_bowtie2_input.index, ch_bowtie2_input.fasta, false, true)
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())

        ch_bam = BOWTIE2_ALIGN.out.bam
        ch_multiqc = ch_multiqc.mix(BOWTIE2_ALIGN.out.log)
    }
    else if (mapper == "bwa") {
        BWA_INDEX(ch_reference)
        ch_versions = ch_versions.mix(BWA_INDEX.out.versions.first())
        ch_bwamem_in = ch_reads
            .join(BWA_INDEX.out.index, by: [0])
            .join(ch_reference, by: [0])
            .multiMap { meta, reads, index, fasta ->
                reads: [meta, reads]
                index: [meta, index]
                fasta: [meta, fasta]
            }


        BWA_MEM(ch_bwamem_in.reads, ch_bwamem_in.index, ch_bwamem_in.fasta, true)
        ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())
        //no mqc for bwa

        ch_bam = BWA_MEM.out.bam
    }
    else {
        error("Unknown mapper: ${mapper}")
    }

    emit:
    bam      = ch_bam       // channel: [ val(meta), [ bam ] ]
    ref      = ch_reference // channel: [ val(meta), [ fasta ] ]
    mqc      = ch_multiqc   // channel: [ multiqc ]
    versions = ch_versions  // channel: [ versions.yml ]
}
