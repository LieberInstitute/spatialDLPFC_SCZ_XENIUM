#!/bin/bash
#SBATCH --job-name=02_get_neuropil_counts
#SBATCH --mem=250G
#SBATCH --time=12:30:00
#SBATCH -n 1
#SBATCH --output=logs/%x_%a.out
#SBATCH --error=logs/%x_%a.out   # file to collect standard output
#SBATCH --array=1-24


echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_CLUSTER_NAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load Modules
module load conda_R/4.5

## List current modules for reproducibility
module list

## Run code
Rscript 02_get_neuropil_counts.R ${SLURM_ARRAY_TASK_ID}


## Memeory stat
#sstat -a -o JobID,MaxVMSizeNode,MaxVMSize,AveVMSize,MaxRSS,AveRS S,MaxDiskRead,MaxDiskWrite,AveCPUFreq,TRESUsageInMax -j ${SLURM_JOB_ID}

echo "**** Job ends ****"
date