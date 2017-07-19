#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=2GB
#SBATCH --job-name=SWM_modelrecov_testmodelfree_leftover
#SBATCH --mail-type=END
#SBATCH --mail-user=aspen.yoo@nyu.edu
#SBATCH --output=modelrecov_testmodelfree_leftover_%a.out

module purge
module load matlab/2016b

cat<<EOF | matlab -nodisplay
addpath(genpath('/home/ay963/matlab-scripts'))
addpath(genpath('/home/ay963/spatialWM'))

expnumber = 2;
blah = num2str($SLURM_ARRAY_TASK_ID);
truemodel = str2double(blah(1));
subjnum = str2double(blah(2));
jobidx = str2double(blah(end-1:end));

nSubj = 11;
testmodel = 2;
runmax = 50;

filepath = ['/home/ay963/spatialWM/fits/exp' num2str(expnumber) '/'];
load([filepath 'modelrecov_truemodel' num2str(truemodel) '_testmodel2_subj' num2str(subjnum+1) '.mat'])

blah = 1:runmax;
blah(runlist_completed) = [];
runlist = blah(jobidx);

fit_parameters(testmodel,subjnum+nSubj+1,runlist,runmax,truemodel,expnumber)

EOF
