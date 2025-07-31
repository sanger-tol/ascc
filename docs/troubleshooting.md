## Troubleshooting

Here we will list some of the common issues we have come across whilst developing ASCC.

## FCS-GX

This error will most likely be solved by the use of `export NXF_SINGULARITY_NEW_PID_NAMESPACE=false` prior to running the nextflow pipeline.

- `exit code 1`
  In cases where the pipeline crashes during the fcsgx process with a very generic `Exited with exit code 1`, please read the `tmp_*` file found in the process working directory. The path should look something like `/PATH/TO/ASCC/work/31/88b821d512d3bd49f74ed684bbe869/tmp_d4305b82-5524-4ef1-9a7e-11061feb5321.summary.txt`

This will not be identical due to the run specific nature of the summary.txt and the work directory so please ensure that you are looking in the right place for your run.

## FILTER_BARCODE

An error of: `ln: failed to create symbolic link` or about not being able to overwrite versions file are typically about the previous process erroring so quickly that the process which manifests the error is a red herring.

In the case of FILTER_BARCODE this is likely due to the use of a non-unique list of barcodes. A simple fix.
