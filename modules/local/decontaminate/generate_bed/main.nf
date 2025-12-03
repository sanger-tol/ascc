process DECONTAMINATE_GENERATE_BED {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81':
        'quay.io/biocontainers/biopython:1.81' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.bed")      , emit: bed
    tuple val(meta), path("*.report")   , emit: report
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = args.ext.prefix   ?: "${meta.id}"
    def args    = args.ext.args     ?: ""
    """
    generate_contamination_bed.py $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        generate_contamination_bed.py: \$(generate_contamination_bed.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix  = args.ext.prefix   ?: "${meta.id}"
    """
    touch ${prefix}_chunked_assembly.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        generate_contamination_bed.py: \$(generate_contamination_bed.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
