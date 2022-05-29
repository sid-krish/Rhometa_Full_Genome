<h1 align="center">Rhometa Full Genome Pipeline</h1>
  <p align="center">
    Reimplementation of LDhat pairwise method for gene conversion as part of developing rhometa

- [About](#about)
- [Getting Started](#getting-started)
  - [Set up using conda](#set-up-using-conda)
  - [Set up using docker](#set-up-using-docker)
- [Quick Start and Output](#quick-start-and-output)
  - [main.nf](#mainnf)
- [Pipeline Options and Advanced Usage](#pipeline-options-and-advanced-usage)
  - [main.nf](#mainnf-1)
- [Issues and Contributing](#issues-and-contributing)
- [Contact](#contact)

## About
Reimplementation of the LDhat pairwise method for gene conversion as part of developing rhometa. For full genome gene conversion type recombination rate analysis it is recommended to use https://github.com/sid-krish/Nextflow_LDhat or the original ldhat pairwise method for gene-conversion https://github.com/auton1/LDhat, however some adjustments need to made get it to working with gene-coversion, this has already been done on Nextflow_LDhat. 

This pipeline will produce identical results when tested with simulated data using https://github.com/sid-krish/Nextflow_LDhat_Sim against Nextflow_LDhat. It was primairly developed to ensure the core algorithms were accurately understood and reimplemented such that it produces exactly the same results when tested with simulated data, however it is not as full featured as the original LDhat pairwise method program and does not have all its functionalities.

The LDhat pairwise method for gene conversion is based on the journal article https://academic.oup.com/genetics/article/160/3/1231/6052507 and implemented in the ldhat package https://github.com/auton1/LDhat. Nextflow_LDhat is simply a pipeline wrapper around LDhat pairwise method configured for gene conversion. Additionaly as in Nextflow_LDhat, the rhometa_full_genome pipeline is also able to perform theta estimates based on equation 1 from the journal article.


<!-- GETTING STARTED -->
## Getting Started

This rhometa_full_genome pipeline is designed to be run on linux and requires nextflow to be installed. 
Dependencies are resolved either via conda or docker images. Support for HPC, docker, singularity, AWS and many other systems are provided via nextflow.

While it is possible to resolve the dependencies using conda for running on macOS, its recommended that this option be used on linux systems for which it has been extensively test.
If running on macOS it recommended that docker be used with the provided image, in which case it is similar to running in a linux environment.

It is also possible to install and run the program on Windows via [wsl](https://docs.microsoft.com/en-us/windows/wsl/install).

### Set up using conda
Instructions for installing nextflow and dependencies via conda
1. Clone the repo
   ```sh
   git clone https://github.com/sid-krish/Rhometa_Full_Genome.git
   ```
2. Install the conda package manager: [Miniconda download](https://conda.io/en/latest/miniconda.html)
3. Install nextflow
   ```sh
   conda install -c bioconda nextflow
   ```
4. Adjust settings in nextflow.config file, by default it is configured to work with docker with modest resources.
   Disable the use of docker by setting the docker option to false. Disabling the use of container engines will cause conda packages to be used by default:
   ```sh
   docker {
       enabled = false
   }
   ```
5. The pipeline is now ready to run, and all dependencies will be automatically resolved with conda.

### Set up using docker
Instructions for installing nextflow and using the provided docker image for dependencies
1. Clone the repo
   ```sh
    git clone https://github.com/sid-krish/Rhometa_Full_Genome.git
   ```
2. Install nextflow [Nextflow install](https://www.nextflow.io/index.html#GetStarted)
3. Install docker desktop [Docker install](https://docs.docker.com/desktop/linux/).
4. Adjust settings in nextflow.config file, by default it is configured to work with docker with modest resources.
5. In the nextflow.config file comment the following line:
   ```
   // conda='environment.yaml'
   ```
6. Ensure docker is running.
7. The pipeline is now ready to run, all the required dependencies are present in the docker image, that the pipeline is preconfigured to use.


<!-- QUICK START AND OUTPUT -->
## Quick Start and Output
The following quick start example makes use of the files in [toy_dataset.zip](https://github.com/sid-krish/Nextflow_LDhat/blob/main/toy_dataset.zip).

This example is designed to help ensure that main.nf (rho estimation pipeline) is configured properly and to demonstrate a typical workflow.
The toy datasets were simulated using the simulation pipeline [Nextflow_LDhat_Sim](https://github.com/sid-krish/Nextflow_LDhat_Sim).

### main.nf
main.nf will by default estimate the theta per site and use it for generating the appropriate lookup table for the number of sequences in the fasta file. The lookup table is then used for the rho estimate.

The default settings of the main.nf pipeline have been configured to work with the quick start example as is.

For this example the command to run is:
```sh
nextflow run main.nf --input_fasta toy_dataset/example.fa
```

The output files will be saved to 'Output'. The file ending with final_results.txt and processed_results.csv have the rho estimates and file ending with theta_est.csv will have the theta (per site) estimate.

The result in the file ending with processed_results.csv should be:
```
max_rho,max_lk
9.0,-41836138.05869978
```

The result in the file ending with theta_est.csv should be:
```
0.010264296769556949
```

<!-- PIPELINE OPTIONS AND ADVANCED USAGE -->
## Pipeline Options and Advanced Usage
### main.nf
In the workflow section of the main.nf script, the following parameters can be adjusted, these are the available options for the pipeline.
```
params.recom_tract_len = 1000
params.lookup_grid = "101,100" // The range of rho values used to generate lookup tables

params.prefix_filename = 'none' // prefix string to output filenames to help distinguish runs
params.input_fasta = 'none'
```

It is also possible to enable the DOWNSAMPLED_LOOKUP_TABLE process to automatically use the downsampled lookup tables generated by lookup_table_gen.nf pipeline in rhometa. To do this the LOOKUP_TABLE_LDPOP process call needs to be commented and the script needs to be adjusted accordingly.


<!-- ISSUES AND CONTRIBUTING -->
## Issues and Contributing
If you have any issues please open an issue with the details and steps for reproducing the issue. If you have any questions please open a issue with the tag "question" or alternatively email one of the authors from the contact section.

If you have a suggestion that would make this better, please fork the repo and create a pull request.


<!-- CONTACT -->
## Contact
Sid Krishnan - sidaswar.krishnan-1@student.uts.edu.au, sid.kr15n@gmail.com \
Aaron Darling - aaron.darling@uts.edu.au \
Matt DeMaere - matthew.demaere@uts.edu.au