process DECONTAMINATE_CLIP_REGIONS_FASTA {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/biopython:1.81'
        : 'quay.io/biocontainers/biopython:1.81'}"

    input:
    tuple val(meta), path(fasta), path(contamination_bed)

    output:
    tuple val(meta), path("*.decontaminated.fa"), emit: decontaminated_fasta
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('clip_regions_bed.py'), eval("clip_regions_bed.py --version | cut -d' ' -f2"), topic: versions, emit: versions_clip_regions_bed

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    clip_regions_bed.py \\
        --sample_id ${prefix} \\
        --fasta ${fasta} \\
        --bed ${contamination_bed}
    """

    stub:
    def prefix = args.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.decontaminated.fa
    """
}
