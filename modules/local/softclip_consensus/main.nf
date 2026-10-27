process SOFTCLIP_CONSENSUS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysamstats:1.1.2--py39he47c912_12'
        : 'quay.io/biocontainers/pysamstats:1.1.2--py311h0152c62_12'}"

    input:
    tuple val(meta), path(fasta, stageAs: 'input/*'), path(bam)

    output:
    tuple val(meta), path("*.fasta")       , emit: fasta
    tuple val(meta), path("*.softclip.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions
    tuple val("${task.process}"), val('pysam'), eval("python -c 'import pysam; print(pysam.__version__)'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    softclip_consensus.py \\
        ${args} \\
        --bam ${bam} \\
        --fasta ${fasta} \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp ${fasta} ${prefix}.fasta
    touch ${prefix}.softclip.tsv
    """
}
