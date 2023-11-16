process FILTER_VECSCREEN_RESULTS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(vecscreen_outfile)

    output:
    tuple val(meta), path("vecscreen.grepped.out"), emit: filtered_vecscreen_outfile
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    VSlistTo1HitPerLine.py --skip_reporting_suspect_hits --skip_reporting_weak_hits --skip_reporting_no_hits ${vecscreen_outfile} > vecscreen.grepped.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        VSlistTo1HitPerLine: \$(VSlistTo1HitPerLine.py -v)
    END_VERSIONS
    """

    stub:
    """
    touch vecscreen.grepped.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        VSlistTo1HitPerLine: \$(VSlistTo1HitPerLine.py -v)
    END_VERSIONS
    """
}
