#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00:00
#SBATCH --mem=2GB
#SBATCH --job-name=SWM_fixed
#SBATCH --mail-type=END
#SBATCH --mail-user=aspen.yoo@nyu.edu
#SBATCH --output=SWM_fixed_%a.out

module purge
module load matlab/2016b

cat<<EOF | matlab -nodisplay
addpath(genpath('/home/ay963/matlab-scripts'))
addpath(genpath('/home/ay963/spatialWM'))

model = 3;
expnumber = 2;
subjnum = $SLURM_ARRAY_TASK_ID;
runlist = 1:50;
runmax = 50;

fit_parameters(model,subjnum,runlist,runmax,[],expnumber)

EOF
