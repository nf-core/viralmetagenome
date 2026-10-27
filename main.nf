#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/viralmetagenome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/viralmetagenome
    Website: https://nf-co.re/viralmetagenome
    Slack  : https://nfcore.slack.com/channels/viralmetagenome
----------------------------------------------------------------------------------------
*/

params.global_prefix = getGlobalPrefix(workflow, params)
def getGlobalPrefix(workflow,params) {
    def date_stamp = new java.util.Date().format( 'yyyyMMdd')
    if (params.prefix) {
        return "${params.prefix}_${date_stamp}_${workflow.manifest.version}_${workflow.runName}".replaceAll("\\s+", "_")
    }
    return null
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VIRALMETAGENOME         } from './workflows/viralmetagenome'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_viralmetagenome_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_viralmetagenome_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_VIRALMETAGENOME {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    VIRALMETAGENOME (
        samplesheet,
        params.outdir,

        // Optional input files
        params.metadata,
        params.blacklist,
        params.contaminants,
        params.adapter_fasta,
        params.spades_yml,
        params.spades_hmm,
        params.mapping_constraints,
        params.annotation_metadata,
        params.multiqc_methods_description,
        params.custom_table_headers,

        // Databases
        params.reference_pool,
        params.kraken2_db,
        params.bracken_db,
        params.kaiju_db,
        params.host_k2_db,
        params.checkv_db,
        params.annotation_db,
        params.prokka_db,

        // Step toggles
        params.skip_preprocessing,
        params.skip_hostremoval,
        params.skip_read_classification,
        params.skip_assembly,
        params.skip_polishing,
        params.skip_precluster,
        params.skip_iterative_refinement,
        params.skip_variant_calling,
        params.skip_vcf_annotation,
        params.skip_consensus_qc,
        params.skip_blast_qc,
        params.skip_checkv,
        params.skip_consensus_annotation,
        params.skip_prokka,

        // Preprocessing
        params.trim_tool,
        params.skip_fastqc,
        params.with_umi,
        params.skip_umi_extract,
        params.umi_discard_read,
        params.umi_deduplicate,
        params.skip_trimming,
        params.save_trimmed_fail,
        params.save_merged,
        params.min_trimmed_reads,
        params.deduplicate,
        params.merge_reads,
        params.skip_complexity_filtering,
        params.decomplexifier,
        params.skip_host_fastqc,
        params.use_host_filtered_reads,

        // Read classification
        params.read_classifiers,
        params.kraken2_save_reads,
        params.kraken2_save_readclassification,
        params.kaiju_taxon_rank,

        // Assembly
        params.assemblers,
        params.normalise_reads,
        params.skip_contig_prinseq,
        params.skip_sspace_basic,
        params.read_distance,
        params.read_distance_sd,
        params.read_orientation,
        params.perc_reads_contig,

        // Contig clustering & polishing
        params.precluster_classifiers,
        params.keep_unclassified,
        params.cluster_method,
        params.identity_threshold,
        params.cluster_with_reference_pool,
        params.skip_singleton_filtering,
        params.min_contig_size,
        params.max_n_perc,
        params.iterative_refinement_cycles,
        params.intermediate_mapper,
        params.call_intermediate_variants,
        params.intermediate_variant_caller,
        params.intermediate_consensus_caller,
        params.intermediate_mapping_stats,

        // Mapping & variant calling
        params.mapper,
        params.variant_caller,
        params.consensus_caller,
        params.mapping_stats,
        params.min_mapped_reads,
        params.keep_unmapped,
        params.ivar_header,

        // Consensus QC
        params.skip_quast,
        params.skip_alignment_qc,

    )
    emit:
    multiqc_report = VIRALMETAGENOME.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden,
        params.merge_reads
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_VIRALMETAGENOME (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_VIRALMETAGENOME.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
