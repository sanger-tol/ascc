process TRAILINGNS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/biopython:1.81'
        : 'quay.io/biocontainers/biopython:1.81'}"

    input:
    tuple val(meta), path(fasta_input_file)

    output:
    tuple val(meta), path("*_trim_Ns"), emit: trailing_ns_report
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('trim_Ns.py'), eval("trim_Ns.py --version | cut -d' ' -f2"), topic: versions, emit: versions_trim_Ns

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = args.ext.prefix ?: "${meta.id}"
    def args = args.ext.args ?: ""
    """
    trim_Ns.py ${fasta_input_file} ${prefix}_trim_Ns ${args}
    """

    stub:
    def prefix = args.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_trim_Ns
    """
}
