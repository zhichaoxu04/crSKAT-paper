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

1. **Job Name**: This assigns a name to the job (or job array).
   1. Single Job
   ```bash
   #BSUB -J job_name
   ```
   2. Job Array: This submits a job array with 123 jobs named "job_array_name".
   ```bash
   #BSUB -J job_array_name[1-123]
   ```
   
   

1. **Output File**:
   1. Output: Redirects the standard output of the job to the specified file.
   ```bash
   #BSUB -o output_file_name
   ```
   2. Error: Redirects the standard error of the job to the specified file.
   ```bash
   #BSUB -e error_file_name
   ```

1. **Requested Resources**:
   - **Number of CPUs**: Requests 40 CPU cores.
     ```bash
     #BSUB -n 40
     ```
   - **Memory**: Requests 128GB of memory.
     ```bash
     #BSUB -M 128
     #BSUB -R rusage[mem=128]
     ```
     
1. **Queue**: Specifies the name of the queue to which the job is submitted (e.g. submit the job to `medium` queue).
   ```bash
   #BSUB -q queue_name (e.g. medium)
   ```
   
1. **Wall Clock Limit**: Specifies the job run time limit in hours and minutes (e.g. 23 hours and 13 minutes).
   ```bash
   #BSUB -W HH:MM (e.g. 23:13)
   ```

1. **Email Notification**: Sends an email when the job starts (`-B` option) and finishes (`-N` option).
   ```bash
   #BSUB -u user_email@example.com
   #BSUB -N 
   #BSUB -B 
   ```

> Note: These are just a selection of the many options available with LSF. Depending on the specific requirements and configuration of your LSF system, you might use a different set of directives. Always refer to your institution's LSF documentation or system administrators for guidance for your own environment.

> Some information was sourced from the official [IBM Spectrum LSF Documentation](https://www.ibm.com/docs/en/spectrum-lsf/10.1.0) available on IBM's website.

## Parallel computing in R script

## Notes:

We acknowledge the support of the High Performance Computing for research facility at the University of Texas MD Anderson Cancer Center ([seadragon](https://fuc.readthedocs.io/en/latest/seadragon.html)) for providing computational resources that have contributed to the research results reported in this paper.
