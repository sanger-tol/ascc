process GENERATE_SAMPLESHEET {
    tag "$meta.id"
    label "process_low"

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'docker.io/ubuntu:20.04' }"

    input:
    tuple val(meta), path(reference)
    path(pacbio_path)
    val(reads_layout)

    output:
    tuple val(meta),    path("samplesheet.csv"),    emit: csv
    path "versions.yml",                            emit: versions

    script:
    def args    = task.ext.args     ?: ""
    def VERSION = "1.1.0"
    """
    echo "Run BTK"
    echo "sample,datatype,datafile,library_layout" > pre_samplesheet.csv

    i=0
    for file in ${pacbio_path}; do
        i=\$((i+1))
        echo "Debug line: Processing file \$file -- T\$i --${reads_layout}"
        echo "${meta.id}_T\$i,pacbio,input_pacbio_files/\$file,${reads_layout}" >> pre_samplesheet.csv
    done

    echo "Debug: MOVE pre_samplesheet.csv to samplesheet.csv"

    mv pre_samplesheet.csv samplesheet.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generate_samplesheet: $VERSION
    END_VERSIONS
    """

    stub:
    def VERSION = "1.1.0"
    """
    touch samplesheet.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generate_samplesheet: $VERSION
    END_VERSIONS
    """
}
