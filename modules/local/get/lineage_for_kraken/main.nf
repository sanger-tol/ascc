process GET_LINEAGE_FOR_KRAKEN {

    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(kraken_file)
    path ncbi_rankedlineage_path

    output:
    path '*_nt_kraken_lineage_file.txt', emit: txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = args.ext.prefix ?: "${meta.id}"
    """
    get_lineage_for_kraken_results.py \\
        $kraken_file \\
        $ncbi_rankedlineage_path \\
        nt \\
        ${prefix}_nt_kraken_lineage_file.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
        general_purpose_functions.py: \$(general_purpose_functions.py --version | cut -d' ' -f2)
        get_lineage_for_kraken_results.py: \$(get_lineage_for_kraken_results.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = args.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_nt_kraken_lineage_file.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(pip list | grep "pandas" | sed 's/[[:blank:]]//g' | sed 's/pandas//g')
        general_purpose_functions.py: \$(general_purpose_functions.py --version | cut -d' ' -f2)
        get_lineage_for_kraken_results.py: \$(get_lineage_for_kraken_results.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
