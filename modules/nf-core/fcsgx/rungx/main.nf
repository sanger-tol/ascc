process FCSGX_RUNGX {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-fcs-gx:0.5.5--h9948957_0':
        'biocontainers/ncbi-fcs-gx:0.5.5--h9948957_0' }"

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
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def module_name = task.ext.module_name ?: ""

    // At Sanger we have a permenant home for the DB on NVME storage
    // def mv_database_to_ram = ramdisk_path ? "rclone copy $gxdb $ramdisk_path/$task.index/" : ''
    // def database = ramdisk_path ? "$ramdisk_path/$task.index/" : gxdb // Use task.index to make memory location unique
    def database = ramdisk_path ?: gxdb

    if ( production_mode ) {
        // Using just the module is not enough
        // Due to how non-user processes set off the module at Sanger
        // We need to create a module config and source it to work
        // see: https://github.com/nextflow-io/nextflow/issues/5980
        """
        echo "Using Production FCSGX with local module"

        modulecmd bash load ${module_name} > .module_def
        source .module_def

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

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fcsgx: \$( gx --help | sed '/build/!d; s/.*:v//; s/;.*//' )
        END_VERSIONS
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

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fcsgx: \$( gx --help | sed '/build/!d; s/.*:v//; s/;.*//' )
        END_VERSIONS
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcsgx: \$( gx --help | sed '/build/!d; s/.*:v//; s/;.*//' )
    END_VERSIONS
    """
}
