process CHUNK_ASSEMBLY_FOR_VECSCREEN {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/biopython:1.81'
        : 'quay.io/biocontainers/biopython:1.81'}"

    input:
    tuple val(meta), path(fasta_input_file)

    output:
    tuple val(meta), path("*_chunked_assembly.fa"), emit: chunked_assembly
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('chunk_assembly_for_vecscreen.py'), eval("chunk_assembly_for_vecscreen.py --version | cut -d' ' -f2"), topic: versions, emit: versions_chunk_assembly_for_vecscreen

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = args.ext.prefix ?: "${meta.id}"
    def args = args.ext.args ?: ""
    """
    chunk_assembly_for_vecscreen.py ${fasta_input_file} ${prefix}_chunked_assembly.fa ${args}
    """

    stub:
    def prefix = args.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_chunked_assembly.fa
    """
}
