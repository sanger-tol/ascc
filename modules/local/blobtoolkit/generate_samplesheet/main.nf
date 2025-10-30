process GENERATE_BTK_SAMPLESHEET {
    tag "$meta.id"
    label "process_low"
    executor 'local'

    input:
    val(meta)
    val(pacbio_path)

    output:
    tuple val(meta),    path("btk_samplesheet.csv"), emit: csv

    exec:
    def fasta_files = files(pacbio_path.resolve('*.fasta.gz'))
    assert fasta_files.size() > 0
    def samplesheet_entries = ["sample,datatype,datafile,library_layout"]
    fasta_files.withIndex().each{ fa, idx -> samplesheet_entries << "${meta.id}_T${idx+1},pacbio,${fa},SINGLE" }
    file(task.workDir.resolve("btk_samplesheet.csv")).text = samplesheet_entries.join("\n")
}
