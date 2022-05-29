<h1 align="center">Rhometa Full Genome Pipeline</h1>
  <p align="center">
    Reimplementation of LDhat pairwise method for gene conversion as part of developing rhometa

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

This example is designed to help ensure that run_ldhat.nf (rho estimation pipeline) and run_theta_est.nf (theta estimation pipeline) are configured properly and to demonstrate a typical workflow.
The toy datasets were simulated using the simulation pipeline [Nextflow_LDhat_Sim](https://github.com/sid-krish/Nextflow_LDhat_Sim).