#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00:00
#SBATCH --mem=2GB
#SBATCH --job-name=SWM_modelrecov_testmodelfixed
#SBATCH --mail-type=END
#SBATCH --mail-user=aspen.yoo@nyu.edu
#SBATCH --output=modelrecov_testmodelfixed_%a.out

module purge
module load matlab/2016b

cat<<EOF | matlab -nodisplay
addpath(genpath('/home/ay963/matlab-scripts'))
addpath(genpath('/home/ay963/spatialWM'))

expnumber = 2;
blah = num2str($SLURM_ARRAY_TASK_ID);
truemodel = str2double(blah(1));
subjnum = str2double(blah(2));
runlistidx = str2double(blah(3));

nSubj = 11;
testmodel = 3;
runmax = 50;
runlist = runlistidx:5:(45+runlistidx);

fit_parameters(testmodel,subjnum+nSubj+1,runlist,runmax,truemodel,expnumber)

EOF
