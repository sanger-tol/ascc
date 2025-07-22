## Troubleshooting

Here we will list some of the common issues we have come across whilst developing ASCC.

## FCS-GX

- `exit code 1`
  In cases where the pipeline crashes during the fcsgx process with a very generic `Exited with exit code 1`, please read the `tmp_*` file found in the process working directory. The path should look something like `/PATH/TO/ASCC/work/31/88b821d512d3bd49f74ed684bbe869/tmp_d4305b82-5524-4ef1-9a7e-11061feb5321.summary.txt`

This will not be identical due to the run specific nature of the summary.txt and the work directory so please ensure that you are looking in the right place for your run.

Towards the bottom of this file you should see the actual error of the run, which could be any of the following.

### Taxonomy not found
