#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:01:00
#SBATCH --mem=1GB
#SBATCH --job-name=test1
#SBATCH --mail-type=END
#SBATCH --mail-user=aspen.yoo@nyu.edu
#SBATCH --output=test2_%j.out

##index=$SLURM_ARRAY_JOB_ID
##job=${SLURM_ARRAY_TASK_ID}
module purge
module load matlab/2016b

cat<<EOF | matlab -nodisplay
addpath(genpath('/home/ay963/matlab-scripts'))
addpath(genpath('/home/ay963/spatialWM'))

expnumber = 2
index = $SLURM_ARRAY_JOB_ID
idx = $SLURM_ARRAY_TASK_ID

slurmtest(expnumber,index,idx)

EOF
