process SELECT_REFERENCE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/16/16fd0599cbc5e52a5ac51f8668ed2c6988b4f44d461606e37953afcd581cd52d/data'
        : 'community.wave.seqera.io/library/biopython_pandas_python:671653bb7f9c4d5b'}"

    input:
    tuple val(meta), path(screen), path(reference), path(reads)

    output:
    tuple val(meta), path("*.json"), path("*_reference.fa"), path(reads), emit: fasta_reads
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    select_reference.py \\
        ${args} \\
        --mash ${screen} \\
        --reference ${reference} \\
        --prefix ${prefix} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(pip show pandas | grep Version | sed 's/Version: //g')
        biopython: \$(pip show biopython | grep Version | sed 's/Version: //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_reference.fa
    touch ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(pip show pandas | grep Version | sed 's/Version: //g')
        biopython: \$(pip show biopython | grep Version | sed 's/Version: //g')
    END_VERSIONS
    """
}
