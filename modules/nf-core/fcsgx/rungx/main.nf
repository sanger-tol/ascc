process FCSGX_RUNGX {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-fcs-gx:0.5.5--h9948957_0':
        'quay.io/biocontainers/ncbi-fcs-gx:0.5.5--h9948957_0' }"

    input:
    tuple val(meta), path(fasta)
    path gxdb
    path ramdisk_path
    val production_mode

    output:
    tuple val(meta), path("*.fcs_gx_report.txt"), emit: fcsgx_report
    tuple val(meta), path("*.taxonomy.rpt")     , emit: taxonomy_report
    tuple val(meta), path("*.summary.txt")      , emit: log
    tuple val(meta), path("*.hits.tsv.gz")      , emit: hits, optional: true
    tuple val("${task.process}"), val('fcsgx'), eval("gx --help | sed '/build/!d; s/.*:v//; s/-.*//'"), emit: versions_fcsgx, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // At Sanger we have a permenant home for the DB on NVME storage
    // def mv_database_to_ram = ramdisk_path ? "rclone copy $gxdb $ramdisk_path/$task.index/" : ''
    // def database = ramdisk_path ? "$ramdisk_path/$task.index/" : gxdb // Use task.index to make memory location unique
    def database = ramdisk_path ?: gxdb

    if ( production_mode ) {
        """
        echo "Using Production FCSGX with local installation"

        export GX_NUM_CORES=${task.cpus}
        export GX_INSTANTIATE_FASTA=1

        run_gx \\
            --fasta ${fasta} \\
            --gx-db ${database} \\
            --tax-id ${meta.taxid} \\
            --generate-logfile true \\
            --out-basename ${prefix} \\
            --out-dir . \\
            ${args}

        """
    } else {
        """
        echo "Using Standard FCSGX with container"

        run_gx.py \\
            --fasta ${fasta} \\
            --gx-db ${database} \\
            --tax-id ${meta.taxid} \\
            --generate-logfile true \\
            --out-basename ${prefix} \\
            --out-dir . \\
            ${args}

        """
    }

    stub:
    // def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fcs_gx_report.txt
    touch ${prefix}.taxonomy.rpt
    touch ${prefix}.summary.txt
    echo "" | gzip > ${prefix}.hits.tsv.gz

    """
}
