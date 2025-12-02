process SOURMASH_SKETCH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sourmash:4.9.4--hdfd78af_0':
        'biocontainers/sourmash:4.9.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(sequence)
    val sketch_params  // Dynamic sketch parameters: "dna -p scaled=X,k=Y,k=Z,..."

    output:
    tuple val(meta), path("*.sig"), emit: signatures
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Use provided sketch_params or fall back to task.ext.args or default
    def args = sketch_params ?: (task.ext.args ?: "dna -p scaled=1000,k=21,k=31,k=51,abund")
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sourmash sketch \\
        $args \\
        --merge '${prefix}' \\
        --output '${prefix}.sig' \\
        $sequence

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix   ?: "${meta.id}"
    """
    touch ${prefix}.sig

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """
}
