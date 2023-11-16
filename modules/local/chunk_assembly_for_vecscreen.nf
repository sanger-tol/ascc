process CHUNK_ASSEMBLY_FOR_VECSCREEN {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::biopython=1.70"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.70--np112py35_1':
        'quay.io/biocontainers/biopython:1.70--np112py35_1' }"

    input:
    tuple val(meta), path(fasta_input_file)

    output:
    tuple val(meta), path("chunked_assembly.fa"), emit: chunked_assembly
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    chunk_assembly_for_vecscreen.py $fasta_input_file chunked_assembly.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        chunk_assembly_for_vecscreen.py: \$(chunk_assembly_for_vecscreen.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    """
    touch chunked_assembly.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        chunk_assembly_for_vecscreen.py: \$(chunk_assembly_for_vecscreen.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
