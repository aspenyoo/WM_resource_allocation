#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --mem=2GB
#SBATCH --job-name=SWM_optimal
#SBATCH --mail-type=END
#SBATCH --mail-user=aspen.yoo@nyu.edu
#SBATCH --output=SWM_optimal_%a.out

module purge
module load matlab/2016b

cat<<EOF | matlab -nodisplay
addpath(genpath('/home/ay963/matlab-scripts'))
addpath(genpath('/home/ay963/spatialWM'))

model = 1;
expnumber = 2;
blah = num2str($SLURM_ARRAY_TASK_ID);
runlistidx = str2double(blah(end-1:end));
subjnum = str2double(blah(1:end-2));

runmax = 50;
runlist = [runlistidx runlistidx+runmax/2];

fit_parameters(model,subjnum,runlist,runmax,[],expnumber)

EOF
