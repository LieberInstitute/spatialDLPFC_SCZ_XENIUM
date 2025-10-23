#!/bin/bash
#SBATCH --job-name=02_find_LR
#SBATCH --mem=100G
#SBATCH --time=10:00:00
#SBATCH -n 1
#SBATCH --output=logs/%x_%a.txt
#SBATCH --error=logs/%x_%a.txt    # file to collect standard output
#SBATCH --array=0-3

ARGS=(neurons neuropil pnn vasc)


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
conda activate ~/mambaforge/envs/scSLAT/
python3 02_find_LR.py --microenv ${ARGS[$SLURM_ARRAY_TASK_ID]}


## Memeory stat
#sstat -a -o JobID,MaxVMSizeNode,MaxVMSize,AveVMSize,MaxRSS,AveRS S,MaxDiskRead,MaxDiskWrite,AveCPUFreq,TRESUsageInMax -j ${SLURM_JOB_ID}

echo "**** Job ends ****"
date