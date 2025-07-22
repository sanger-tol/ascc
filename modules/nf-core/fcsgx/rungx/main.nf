process FCSGX_RUNGX {
    tag "$meta.id"
    label 'process_high'
    maxForks 1          // Should ensure that only one instance runs at any given time
    maxRetries 2        // 2 retries to hopefully avoid the issue of fcs crashes once and kills the pipeline
    errorStrategy { sleep(1200 as long); return 'retry' } // 20 Minute delay, same reason as above.


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-fcs-gx:0.5.5--h9948957_0':
        'biocontainers/ncbi-fcs-gx:0.5.5--h9948957_0' }"

    input:
    tuple val(meta), path(fasta)
    path gxdb
    path ramdisk_path

    output:
    tuple val(meta), path("*.fcs_gx_report.txt"), emit: fcsgx_report
    tuple val(meta), path("*.taxonomy.rpt")     , emit: taxonomy_report
    tuple val(meta), path("*.summary.txt")      , emit: log
    tuple val(meta), path("*.hits.tsv.gz")      , emit: hits, optional: true
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // At Sanger we have a permenant home for the DB on NVME storage
    // def mv_database_to_ram = ramdisk_path ? "rclone copy $gxdb $ramdisk_path/$task.index/" : ''
    // def database = ramdisk_path ? "$ramdisk_path/$task.index/" : gxdb // Use task.index to make memory location unique
    def database = ramdisk_path ?: gxdb
    // cp `readlink` has been used as a potentially better alternative
    // to to adding the below to the script block.
    // In theory they should do the same thing.
    //    export GX_INSTANTIATE_FASTA=1

    """
    export GX_NUM_CORES=${task.cpus}
    export FCS_DEFAULT_IMAGE="https://depot.galaxyproject.org/singularity/ncbi-fcs-gx:0.5.5--h9948957_0"

    cp `readlink ${fasta}` ${prefix}_copied_input.fasta

    run_gx.py \\
        --fasta ${prefix}_copied_input.fasta \\
        --gx-db ${database} \\
        --tax-id ${meta.taxid} \\
        --generate-logfile true \\
        --out-basename ${prefix} \\
        --out-dir . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcsgx: \$( gx --help | sed '/build/!d; s/.*:v//; s/;.*//' )
    END_VERSIONS
    """

    stub:
    // def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fcs_gx_report.txt
    touch ${prefix}.taxonomy.rpt
    touch ${prefix}.summary.txt
    echo "" | gzip > ${prefix}.hits.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcsgx: \$( gx --help | sed '/build/!d; s/.*:v//; s/;.*//' )
    END_VERSIONS
    """
}
