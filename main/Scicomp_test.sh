#!/bin/bash -l
#$ -pe smp 2
#$ -cwd 
#$ -o log.out
#$ -e log.err
#$ -q all.q

ml nextflow/24.04.2

nextflow run scicomp_main.nf \
  -c scicomp_fullWorkflow.config \
  --publish_dir Output/ \
  --source ncbi \
  --reads input.csv 
