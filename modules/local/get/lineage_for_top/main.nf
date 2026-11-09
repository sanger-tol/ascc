process GET_LINEAGE_FOR_TOP {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9'
        : 'biocontainers/python:3.9'}"

    input:
    tuple val(meta), path(tophits)
    path ncbi_lineage_path

    output:
    tuple val(meta), path("*.csv"), emit: full
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('get_lineage_for_top.py'), eval("get_lineage_for_top.py --version | sed 's/^//'"), topic: versions, emit: versions_get_lineage_for_top

    script:
    """
    get_lineage_for_top.py --blast_tsv ${tophits} --taxdump ${ncbi_lineage_path} --output_csv ./${meta.id}_BLAST_results_with_lineage.csv --column_name_prefix nt_
    """

    stub:
    """
    touch ${meta.id}_BLAST_results_with_lineage.csv
    """
}
