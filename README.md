## crSKAT-paper

[crSKAT](https://github.com/zhichaoxu04/crSKAT): A powerful tests for genetic set-based inference with interval-censored competing risks outcomes.

We provide the scripts to perform:

- Type I error simulations

- Power simulations

- Real data analysis with [UK Biobank data](https://www.ukbiobank.ac.uk)

- Diagnostic plots


## LSF files
We provide the LSF files to submit R jobs for simulations and real data analysis using parallel computing. 
When writing an LSF file (often a batch script), you'll use the `#BSUB` directive to specify job parameters.

1. **Job Name**: This assigns a name to the job.
   ```bash
   #BSUB -J job_name(e.g. T1Error)
2. **Output Files**:
   1. Output: Redirects the standard output of the job to the specified file.
   ```bash
   #BSUB -o output_file_name (e.g. /your/path/to/Log.out)
   ```
   2. Error: Redirects the standard error of the job to the specified file.
   ```bash
   #BSUB -o output_file_name (e.g. /your/path/to/Error.err)
   ```
4. **Output File**:
   ```bash
   #BSUB -o output_file_name (e.g. /your/path/to/Log.out)
5. **Output File**:
   ```bash
   #BSUB -o output_file_name (e.g. /your/path/to/Log.out)
6. **Output File**:
   ```bash
   #BSUB -o output_file_name (e.g. /your/path/to/Log.out)
7. **Output File**:
   ```bash
   #BSUB -o output_file_name (e.g. /your/path/to/Log.out)

## Parallel computing in R script

## Notes:

We acknowledge the support of the High Performance Computing for research facility at the University of Texas MD Anderson Cancer Center ([seadragon](https://fuc.readthedocs.io/en/latest/seadragon.html)) for providing computational resources that have contributed to the research results reported in this paper.
