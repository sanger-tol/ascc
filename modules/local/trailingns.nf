process TRAILINGNS {
    tag "$meta.id"
    label 'process_single'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.70--np112py27_1':
        'biocontainers/biopython:1.70--np112py27_1' }"

    input:
    tuple val(meta), path(fasta_input_file)

    output:
    tuple val(meta), path("*_trim_Ns"), emit: trailing_ns_report
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = args.ext.prefix   ?: "${meta.id}"
    def args    = args.ext.args     ?: ""
    """
    trim_Ns.py $fasta_input_file ${prefix}_trim_Ns ${args}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        trim_Ns.py: \$(trim_Ns.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix  = args.ext.prefix   ?: "${meta.id}"
    def args    = args.ext.args     ?: ""
    """
    touch ${prefix}_trim_Ns
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        trim_Ns.py: \$(trim_Ns.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}