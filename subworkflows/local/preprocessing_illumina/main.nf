// modules

include { lowReadSamplesToMultiQC            } from '../utils_nfcore_viralmetagenome_pipeline'
include { PRINSEQPLUSPLUS as PRINSEQ_READS   } from '../../../modules/nf-core/prinseqplusplus/main'
include { HUMID                              } from '../../../modules/nf-core/humid/main'
include { BBMAP_BBDUK                        } from '../../../modules/nf-core/bbmap/bbduk/main'
include { CAT_FASTQ                          } from '../../../modules/nf-core/cat/fastq/main'
include { FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC  } from '../fastq_fastqc_umitools_trimmomatic'
include { FASTQ_FASTQC_UMITOOLS_FASTP        } from '../../nf-core/fastq_fastqc_umitools_fastp/main'
include { FASTQ_KRAKEN_HOST_REMOVE           } from '../fastq_kraken_host_remove'

workflow PREPROCESSING_ILLUMINA {

    take:
    ch_reads                   // channel: [ [ meta ], [ ch_reads ] ]
    ch_kraken2_host_db         // channel: [ path(kraken2_host_db) ]
    ch_adapter_fasta           // channel: [ path(adapter_fasta) ]
    ch_contaminants            // channel: [ path(contaminants_fasta) ]
    trim_tool                  // string:  [ fastp | trimmomatic ]
    skip_fastqc                // boolean: skip FastQC on raw and trimmed reads
    with_umi                   // boolean: reads carry UMIs
    skip_umi_extract           // boolean: skip UMI extraction
    umi_discard_read           // integer: [ 0 | 1 | 2 ] read to discard after UMI extraction
    skip_trimming              // boolean: skip adapter/quality trimming
    save_trimmed_fail          // boolean: fastp, save reads failing trimming
    save_merged                // boolean: fastp, save merged read pairs
    min_trimmed_reads          // integer: minimum reads after trimming to keep a sample
    umi_deduplicate            // string:  [ read | mapping | both ] where UMI deduplication happens
    deduplicate                // boolean: deduplicate reads
    merge_reads                // boolean: concatenate reads belonging to the same sample / group
    skip_complexity_filtering  // boolean: skip low-complexity read filtering
    decomplexifier             // string:  [ bbduk | prinseq ]
    skip_hostremoval           // boolean: skip host read removal
    skip_host_fastqc           // boolean: skip FastQC after host removal

    main:
    ch_multiqc_files    = channel.empty()
    ch_trim_read_count  = channel.empty()

    // QC & UMI & Trimming with fastp or trimmomatic
    if (trim_tool == 'trimmomatic') {
        FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC (
            ch_reads,
            skip_fastqc,
            with_umi,
            skip_umi_extract,
            umi_discard_read,
            skip_trimming,
            min_trimmed_reads
            )
        ch_trim_read_count = ch_trim_read_count.mix(FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC.out.trim_read_count)
        ch_multiqc_files   = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC.out.fastqc_raw_zip)
        ch_multiqc_files   = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC.out.fastqc_trim_html)
        ch_multiqc_files   = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC.out.trim_log)
        ch_multiqc_files   = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC.out.umi_log)

        ch_reads_trim = FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC.out.reads
    }
    else if (trim_tool == 'fastp') {
        fastp_reads = ch_reads.map{meta, reads -> [meta, reads, ch_adapter_fasta] }
        FASTQ_FASTQC_UMITOOLS_FASTP (
            fastp_reads,
            skip_fastqc,
            with_umi,
            skip_umi_extract,
            umi_discard_read,
            skip_trimming,
            save_trimmed_fail,
            save_merged,
            min_trimmed_reads
            )

        ch_trim_read_count = ch_trim_read_count.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_read_count)
        ch_multiqc_files   = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_raw_zip)
        ch_multiqc_files   = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_trim_zip)
        ch_multiqc_files   = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_json)
        ch_multiqc_files   = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.umi_log)


        ch_reads_trim = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads
    }

    // Keeping track of failed reads for reporting
    ch_failed_reads = ch_trim_read_count
        .filter{_meta, num_reads -> num_reads < min_trimmed_reads.toLong() }

    // deduplicate UMI's with HUMID
    if (with_umi && ['read', 'both'].contains(umi_deduplicate) && deduplicate ) {
        HUMID (
            ch_reads_trim,
            [[:],[]]
        )
        ch_reads_dedup   = HUMID.out.dedup
        ch_multiqc_files = ch_multiqc_files.mix(HUMID.out.stats)
    }
    else {
        ch_reads_dedup = ch_reads_trim
    }

    // Merge reads belonging to the same sample / group
    if (merge_reads) {
        ch_reads_grouped = ch_reads_dedup
            .map { meta, reads -> [meta + [id: meta.sample], reads] }
            .groupTuple ()

        CAT_FASTQ ( ch_reads_grouped.map { meta, reads -> [meta, reads.flatten()] } )
        ch_reads_dedup_joined = CAT_FASTQ.out.reads
    } else {
        ch_reads_dedup_joined = ch_reads_dedup
    }


    // Decomplexification with BBDuk
    if (!skip_complexity_filtering) {
        if (decomplexifier == 'bbduk') {
            BBMAP_BBDUK (
                ch_reads_dedup_joined,
                ch_contaminants,
            )
            ch_reads_decomplexified = BBMAP_BBDUK.out.reads
            ch_multiqc_files        = ch_multiqc_files.mix(BBMAP_BBDUK.out.log)
        } else if (decomplexifier == 'prinseq') {
            ch_prinseq_in = ch_reads_dedup_joined.map { meta, reads -> [meta, reads, []] }
            PRINSEQ_READS (
                ch_prinseq_in
            )
            ch_reads_decomplexified = PRINSEQ_READS.out.good_reads
            ch_multiqc_files        = ch_multiqc_files.mix(PRINSEQ_READS.out.log)
        }
    } else {
        ch_reads_decomplexified = ch_reads_dedup_joined
    }

    // Host removal with kraken2
    if (!skip_hostremoval){
        FASTQ_KRAKEN_HOST_REMOVE (
            ch_reads_decomplexified,
            ch_kraken2_host_db,
            skip_host_fastqc,
            min_trimmed_reads,
        )

        ch_reads_hostremoved   = FASTQ_KRAKEN_HOST_REMOVE.out.reads_hostremoved
        ch_failed_reads        = ch_failed_reads.mix(FASTQ_KRAKEN_HOST_REMOVE.out.reads_hostremoved_fail)
        ch_multiqc_files       = ch_multiqc_files.mix( FASTQ_KRAKEN_HOST_REMOVE.out.mqc )

    } else {
        ch_reads_hostremoved = ch_reads_decomplexified
    }


    //
    // Create a section that reports failed samples and their read counts
    //
    ch_low_reads_mqc = lowReadSamplesToMultiQC(ch_failed_reads, min_trimmed_reads)
        .collectFile(name:'samples_low_reads_mqc.tsv')


    emit:
    reads                   = ch_reads_hostremoved            // channel: [ [ meta ], [ reads ] ]
    reads_decomplexified    = ch_reads_decomplexified         // channel: [ [ meta ], [ reads ] ]
    reads_trimmed           = ch_reads_dedup                  // channel: [ [ meta ], [ reads ] ]
    mqc                     = ch_multiqc_files                // channel: [ [ meta ], [ mqc ] ]
    low_reads_mqc           = ch_low_reads_mqc                // channel: [ mqc ]
}
