## crSKAT-paper

[crSKAT](https://github.com/zhichaoxu04/crSKAT): A powerful tests for genetic set-based inference with interval-censored competing risks outcomes.

We provide the scripts to perform:

- Type I error simulations

- Power simulations

- Real data analysis with [UK Biobank data](https://www.ukbiobank.ac.uk)

- Diagnostic plots

> [!IMPORTANT]  
> For more details or instructions of this package, please refer to [crSKAT](https://github.com/zhichaoxu04/crSKAT).

## Getting Started

Download and install following required R packages:

- Download [crSKAT](https://github.com/zhichaoxu04/crSKAT) package from
  Github using:

<!-- -->

    git clone https://github.com/zhichaoxu04/crSKAT.git

- Or, install [crSKAT](https://github.com/zhichaoxu04/crSKAT) package in
  R directly

  - First, install [devtools](https://devtools.r-lib.org) in R from
    CRAN:

    ``` r
    install.packages("devtools")
    ```

  - Then, install [crSKAT](https://github.com/zhichaoxu04/crSKAT) using
    the `install_github` function and load the package:

    ``` r
    devtools::install_github("zhichaoxu04/crSKAT")
    library(crSKAT)
    ```

- Make sure that all the required packages have been installed or
  updated. Here are some of the required packages:

  - [CompQuadForm](https://cran.r-project.org/web/packages/CompQuadForm/index.html):
    Compute the distribution function of quadratic forms in normal
    variables using Imhof’s method, Davies’s algorithm, Farebrother’s
    algorithm or Liu et al.’s algorithm.
  - [nleqslv](https://cran.r-project.org/web/packages/nleqslv/index.html):
    Solves a system of nonlinear equations using a Broyden or a Newton
    method with a choice of global strategies such as line search and
    trust region.
  - [ICSKAT](https://cran.r-project.org/web/packages/ICSKAT/index.html):
    Implements the Interval-Censored Sequence Kernel Association
    (ICSKAT) test for testing the association between interval-censored
    time-to-event outcomes and groups of single nucleotide polymorphisms
    (SNPs).
  - [bindata](https://cran.r-project.org/web/packages/bindata/index.html):
    Generates correlated artificial binary data.

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
Parallel computing in R allows you to take advantage of multiple cores or processors to speed up computations. There are several packages and techniques available in R to facilitate parallel processing. Here, I'll provide an example using the `foreach` and `doParallel` R packages in our paper:

## Notes:

We acknowledge the support of the High Performance Computing for research facility at the University of Texas MD Anderson Cancer Center ([seadragon](https://fuc.readthedocs.io/en/latest/seadragon.html)) for providing computational resources that have contributed to the research results reported in this paper.
