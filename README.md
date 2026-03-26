# Pulsenet 2.0 Bioinformatics Pipeline

## Introduction

This repo contains the source code for the bioinformatics pipeline used in the Pulsenet 2.0 application. The pipeline is organized into two directories, "main/" and "processes/", which contain the source code for the main nextflow pipeline (main/) and the code and container recipes used for each process in that pipeline (processes/).

The code in this repo is provided as is and with no warranty or support.

## Using this repo

Briefly, to use this repo you would need to build the process containers used by the main pipeline, set up your nextflow configuration files, and then run one of the nextflow entrypoint scripts.

### pipeline code in the main/ directory

Within the main/ directory are the following important contents:

1. two entrypoint scripts, main.nf and scicomp_main.nf. These are entrypoints used for our cloud and HPC environments, respectively.
2. scicomp_fullWorkflow.config, which is the nextflow configuration file used to run scicomp_main.nf in our HPC environment.
3. a directory named modules/ which contains nextflow scripts for the processes used in the PN2.0 pipeline.
4. Other nextflow code and database files used in the pipeline.

As specified in the scicomp_fullWorkflow.config file, each process uses a container for its execution. The recipes to build each of these containers is provided in the processes/ directory. Each nextflow script in the main/workflows/ directory corresponds to a subdirectory within processes/ with the same name. For example, the process file AlleleCalling.nf is used with a container that can be built from the files within the directory as described below.

### Container recipes in the processes/ directory

Within each subdirectory in processes/ is a Dockerfile and any other required files to build a container for the module script with the same name as the directory. Any nextflow modules that do not have a corresponding processes/ subdirectory use publicly available containers.

