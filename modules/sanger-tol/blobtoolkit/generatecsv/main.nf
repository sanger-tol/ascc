process BLOBTOOLKIT_GENERATECSV {
    tag "${meta.id ?: meta2.id ?: meta3.id}"
    executor "local"

    input:
    tuple val(meta),    path(pacbio_files)
    tuple val(meta2),   path(ont_files)
    tuple val(meta3),   path(illumina_files), val(illumina_layout)

    output:
    tuple val(new_meta),    path("${prefix}.samplesheet.csv")   , emit: csv
    path("versions.yml")                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    exec:
    // Note: Manually bump version number when updating module
    def VERSION       = "1.0.0"

    // CHECK META.ID, IF 1 NOT SET THEN TRY NEXT
    sample_id = meta.id ?: meta2.id ?: meta3.id
    new_meta = [id: sample_id]
    prefix = task.ext.prefix ?: sample_id

    samplesheet_entries =
      [["sample", "datatype", "datafile", "library_layout"].join(",")] +
      [
        [meta  + [data_type: "pacbio"], pacbio_files,   "SINGLE"],
        [meta2 + [data_type: "ont"],    ont_files,      "SINGLE"],
        [meta3 + [data_type: "hic"],    illumina_files, illumina_layout]
      ].collect {
        this_meta, lst, layout ->
        lst.collect {
          fa ->
          [this_meta.data_type, fa, layout].join(",")
        }
      }
      .flatten()
      .withIndex()
      .collect { str, idx ->
        "${prefix}_T${idx+1}," + str
      }

    file(task.workDir.resolve("${prefix}.samplesheet.csv")).text = samplesheet_entries.join("\n")

    file("${task.workDir}/versions.yml").text = """\
        BLOBTOOLKIT_GENERATECSV:
            blobtoolkit_generatecsv: ${VERSION}
        """.stripIndent()
}
