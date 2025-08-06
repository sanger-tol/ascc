## Troubleshooting

Here we will list some of the common issues we have come across whilst developing ASCC.

## FCS-GX

This error will most likely be solved by the use of `export NXF_SINGULARITY_NEW_PID_NAMESPACE=false` prior to running the nextflow pipeline.

- `exit code 1`
  In cases where the pipeline crashes during the fcsgx process with a very generic `Exited with exit code 1`, please read the `tmp_*` file found in the process working directory. The path should look something like `/PATH/TO/ASCC/work/31/88b821d512d3bd49f74ed684bbe869/tmp_d4305b82-5524-4ef1-9a7e-11061feb5321.summary.txt`

This will not be identical due to the run specific nature of the summary.txt and the work directory so please ensure that you are looking in the right place for your run.

## FCS-GX part 2

If the above fix does not work and FCS continues to crash in what seems like a generic fashion, try running fcs locally. With NO container.

To make life easier, there is now `--fcs_override_samplesheet` (A boolean, default false) and `--fcs_samplesheet` (For the samplesheet path). The samplesheet can be generated with the script `./bin/ascc_fcsgx_wrapper.py`. Point this towards a local install of fcs along with a number of other config options that you can read about using `-h`.

This will parse the ASCC samplesheet, gunzip input if needed, filter the fasta (using the same method as in ASCC), run fcs-gx, parse the fcs-gx results into something ASCC can read and fianlly generate the fcs_samplesheet.csv file needed for the `--fcs_override_samplesheet` flag.
