import java.nio.file.Paths
import java.nio.file.Files

process NEXTFLOW_RUN {
    tag "$pipeline_name"

    input:
    val pipeline_name     // String
    val nextflow_opts     // String
    val nextflow_files    // Map [ params-file: params.yml , c: configs/multiqc.config ]
    val pipeline_files    // Map [ input: samplesheet.csv ]
    val additional_config // custom configs

    when:
    task.ext.when == null || task.ext.when

    exec:
    // def args = task.ext.args ?: ''
    def cache_dir = Paths.get(workflow.workDir.resolve(pipeline_name).toUri())
    Files.createDirectories(cache_dir)
    def nxf_cmd = [
        'nextflow run',
            pipeline_name,
            nextflow_opts,
            nextflow_files ? nextflow_files.collect{ key, value -> "-$key $value" }.join(' ') : '',
            pipeline_files ? pipeline_files.collect{ key, value -> "--$key $value" }.join(' ') : '',
            "--outdir $task.workDir/results",
    ]

    ProcessBuilder builder = new ProcessBuilder(nxf_cmd.join(" ").tokenize(" "))
    builder.directory(cache_dir.toFile())
    def process = builder.start()

    // Read stdout and stderr concurrently
    def output_data = new StringBuilder()
    def error = new StringBuilder()

    def stdoutThread = Thread.start {
        process.inputStream.eachLine { line -> output_data.append(line).append("\n") }
    }
    def stderrThread = Thread.start {
        process.errorStream.eachLine { line -> error.append(line).append("\n") }
    }

    // Wait for the process to complete and join threads
    def exitCode = process.waitFor()

    stdoutThread.join()
    stderrThread.join()

    // Check the exit code
    assert exitCode == 0 : "Pipeline failed with exit code ${exitCode}\nError: ${error}\nOutput: ${output_data}"

    // Emit results
    output:
    path "results", emit: output
    //val output_data, emit: log // <-- This need investigating, why is the output_data not assigned at this point but is in the original version by mahesh?
}
