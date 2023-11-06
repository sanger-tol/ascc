process CUSTOM_GETCHROMSIZES {
    tag "$fasta"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
    tuple   val(meta), path(fasta, stageAs: 'input/*')
    val     suffix

    output:
    tuple val(meta), path ("*.${suffix}")   , emit: sizes
    tuple val(meta), path ("*.fa")          , emit: fasta
    tuple val(meta), path ("*.fai")         , emit: fai
    tuple val(meta), path ("*.gzi")         , emit: gzi, optional: true
    path  "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """

    ln -s ${fasta} ${prefix}.fa
    samtools faidx ${prefix}.fa -o ${prefix}.fa.fai
    cut -f 1,2 ${prefix}.fa.fai > ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getchromsizes: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${prefix}.fa
    touch ${prefix}.fa.fai
    touch ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getchromsizes: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
