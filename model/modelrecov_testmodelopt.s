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
#SBATCH --output=modelrecov_testmodelopt_%a.out

module purge
module load matlab/2016b

cat<<EOF | matlab -nodisplay
addpath(genpath('/home/ay963/matlab-scripts'))
addpath(genpath('/home/ay963/spatialWM'))

expnumber = 2;
blah = num2str($SLURM_ARRAY_TASK_ID);
truemodel = str2double(blah(1));
subjnum = str2double(blah(2:end-2));
runlistidx = str2double(blah(end-1:end));

testmodel = 1;
runmax = 50;
runlist = [runlistidx runlistidx+runmax/2];

fit_parameters(model,subjnum,runlist,runmax,[],expnumber)

EOF
