#!/bin/bash
#SBATCH --job-name=06_plot_clusters_on_tissue
#SBATCH --mem=200G
#SBATCH --time=20:00:00
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
Rscript 06_plot_clusters_on_tissue.R


## Memeory stat
#sstat -a -o JobID,MaxVMSizeNode,MaxVMSize,AveVMSize,MaxRSS,AveRS S,MaxDiskRead,MaxDiskWrite,AveCPUFreq,TRESUsageInMax -j ${SLURM_JOB_ID}

echo "**** Job ends ****"
date