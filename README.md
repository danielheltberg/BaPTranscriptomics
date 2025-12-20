# BaPTranscriptomics
Script used in 'Synergetic and independent mechanisms of generating plasticity in response to benzo[a]pyrene exposure in Caenorhabditis elegans and Gadus morhua'
# Scripts are either in Bash for linux command line processing of raw RNA reads OR in R scripts/txt files. Both will have modifications for either C elegans or G morhua
You may wish to submit jobs to a workload manager such as slurm. We used the following code to do so but omitted it in the workflow
*#!/bin/bash*
*#SBATCH -c 8 --mem 64G*
