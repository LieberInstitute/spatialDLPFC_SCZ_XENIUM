#!/bin/bash
#SBATCH --job-name=06_layer_stratified_de
#SBATCH --mem=200G
#SBATCH --time=30:00
#SBATCH -n 1
#SBATCH --output=logs/%x.txt
#SBATCH --error=logs/%x.txt    # file to collect standard output


echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_CLUSTER_NAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load Modules
module load conda_R/4.3.x

## List current modules for reproducibility
module list

## Run code
Rscript 06_layer_stratified_de.R

echo "**** Job ends ****"
date