process DECONTAMINATE_CLIP_REGIONS_FASTA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81':
        'quay.io/biocontainers/biopython:1.81' }"

    input:
    tuple val(meta), path(fasta_input_file)
    path(contamination_bed)

    output:
    tuple val(meta), path("*.decontaminated.fasta")  , emit: decontaminated_fasta
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = args.ext.prefix   ?: "${meta.id}"
    def args    = args.ext.args     ?: ""
    """
    clip_regions_bed.py $fasta_input_file $contamination_bed ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        clip_regions_bed.py: \$(clip_regions_bed.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix  = args.ext.prefix   ?: "${meta.id}"
    """
    touch ${prefix}.decontaminated.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        clip_regions_bed.py: \$(clip_regions_bed.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
