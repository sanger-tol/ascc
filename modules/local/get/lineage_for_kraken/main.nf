process GET_LINEAGE_FOR_KRAKEN {

    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pandas:1.5.2'
        : 'quay.io/biocontainers/pandas:1.5.2'}"

    input:
    tuple val(meta), path(kraken_file)
    path ncbi_rankedlineage_path

    output:
    tuple val(meta), path('*_nt_kraken_lineage_file.txt'), emit: txt
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pandas'), eval("python3 -c 'import pandas; print(pandas.__version__)'"), topic: versions, emit: versions_pandas
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('general_purpose_functions.py'), val("1.0.0"), topic: versions, emit: versions_general_purpose_functions
    tuple val("${task.process}"), val('get_lineage_for_kraken_results.py'), eval("get_lineage_for_kraken_results.py --version | cut -d' ' -f2"), topic: versions, emit: versions_get_lineage_for_kraken_results

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = args.ext.prefix ?: "${meta.id}"
    """
    get_lineage_for_kraken_results.py \
        ${kraken_file} \
        ${ncbi_rankedlineage_path} \
        nt \
        ${prefix}_nt_kraken_lineage_file.txt
    """

    stub:
    def prefix = args.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_nt_kraken_lineage_file.txt
    """
}
