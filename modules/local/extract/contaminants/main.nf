process EXTRACT_CONTAMINANTS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-7fa7184beddbf01b0b0ae64ed643f6c05f12cbcc:337918d8410b938b17fd4beb45024a78ffa6b0d3-0'
        : 'biocontainers/mulled-v2-7fa7184beddbf01b0b0ae64ed643f6c05f12cbcc:337918d8410b938b17fd4beb45024a78ffa6b0d3-0'}"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(blast_data)

    output:
    tuple val(meta), path("*.ALL*unfiltered_scaffold_coverage.bed"), emit: contamination_bed
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('extract_contaminants_by_type.py'), eval("extract_contaminants_by_type.py -v"), topic: versions, emit: versions_extract_contaminants_by_type

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    extract_contaminants_by_type.py \
        ${blast_data} \
        --assembly_file ${fasta} \
        --out_prefix ${prefix}
    """

    stub:
    """
    touch test.bed
    """
}
