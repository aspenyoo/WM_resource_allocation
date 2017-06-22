#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --mem=2GB
#SBATCH --job-name=SWM_optimal
#SBATCH --mail-type=END
#SBATCH --mail-user=aspen.yoo@nyu.edu
#SBATCH --output=SWM_optimal_%j.out

index=${PBS_ARRAYID}
job=${PBS_JOBID}
walltime_lim=${PBS_WALLTIME}
script_name=${PBS_JOBNAME}
module purge
module load matlab/2016b

export MATLABPATH=/home/ay963/matlab-scripts
cat<<EOF | matlab -nodisplay
addpath(genpath('/home/ay963/matlab-scripts'))
addpath(genpath('/home/ay963/spatialWM'))

expnumber = 2;
blah = num2str($index);
runlistidx = str2double(blah(1));
subjnum = str2double(blah(2:end));

rng(runlistidx)
runmax = 50;
runlist = runlistidx:5:(45+runlistidx);
runlist(randperm(10)) = runlist;

fit_parameters(1,subjnum,runlist,runmax,[],expnumber)

EOF
