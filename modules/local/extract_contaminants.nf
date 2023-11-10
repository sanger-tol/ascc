process EXTRACT_CONTAMINANTS {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9 conda-forge::biopython=1.78 conda-forge::pybedtools=0.9.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'biocontainers/biopython:1.78' }"

    input:
    tuple val(meta), path(blast_data)
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path( "assembly.ALL.unfiltered_scaffold_coverage.bed" ), emit: contamination_bed
    path "versions.yml"                                                     , emit: versions

    script:
    def args    = task.ext.args ?: ''
    """
    extract_contaminants_by_type.py \\
        ${blast_data} \\
        --assembly_file ${fasta} \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        extract_contaminants_by_type: \$(extract_contaminants_by_type.py -v)
    END_VERSIONS
    """

    stub:
    def args    = task.ext.args ?: ''
    """
    touch test.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        extract_contaminants_by_type: \$(extract_contaminants_by_type.py -v)
    END_VERSIONS
    """
}
