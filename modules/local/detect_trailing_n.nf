process DETECT_TRAILING_N {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::biopython=1.70"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.70--np112py35_1':
        'quay.io/biocontainers/biopython:1.70--np112py35_1' }"

    input:
    tuple val(meta), path(fasta_input_file)

    output:
    tuple val(meta), path("*_trim_Ns"), emit: trailing_n_report
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    trim_Ns.py ${fasta_input_file} ${prefix}_trim_Ns
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        trim_Ns.py: \$(trim_Ns.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_trim_Ns

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        trim_Ns.py: \$(trim_Ns.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}