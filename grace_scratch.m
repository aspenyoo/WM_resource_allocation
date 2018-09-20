%% load data for grace's experiment

clear all; clc

load('data/grace_data/n_9_all_data.mat');
% load('data/grace_data/raw/noTMS/cleandata.mat')

exppriorityVec = [2/3 1/3];


%% GET DATA IN CORRECT FORMAT FOR FITTING

% error
error_i = all_data.s_all.i_sacc_err;
error_f = all_data.s_all.f_sacc_err;

% indicate which subject, condition, and hemifield to use
subjnum = 2; % 1-9
TMScond = 1; % 1 = no tms, 2 = ips2, 3 = spcs
hemifield = 2; % 1: left hemi. 2: right hemi. 0: both hemi

% get data from a particular subject and condition (and hemifield?)
idx = all_data.use_trial;
idx = idx & (all_data.subj_all == subjnum);
idx = idx & (all_data.TMS_cond_all == TMScond);
if (hemifield); idx = idx & (all_data.s_all.trialinfo(:,2) == hemifield); end

data = cell(1,2);
idxx = idx & (all_data.s_all.trialinfo(:,1) == 31); % high priority
data{1} = error_i(idxx,:);
idxx = idx & (all_data.s_all.trialinfo(:,1) == 32); % low priority
data{2} = error_f(idxx,:);

% delete any nans
data{1}(isnan(data{1})) = [];
data{2}(isnan(data{2})) = [];

data

%% FIT MODEL

filepath = 'fits/grace/';
model = 'flexible';
runlist = 2:50;
runmax = 50;
fixparams = [];

filename = [filepath 'fits_model_' num2str(model) '_subj' num2str(subjnum) ...
    '_TMScond_' num2str(TMScond) 'hemi' num2str(hemifield) '.mat'];

try load(filename); catch; ML_parameters = []; nLLVec = []; end
try runlist_completed*2; catch; runlist_completed = []; end % seeing if runlist_completed exists yet
   

for irun = 1:length(runlist);
    runlistt = runlist(irun);
    try
        [bfp,fval,rlc] = fit_parameters(model,data,exppriorityVec,runlistt,runmax,fixparams);
        
        ML_parameters = [ML_parameters; bfp];
        nLLVec = [nLLVec fval];
        runlist_completed = [runlist_completed rlc];
        save(filename,'ML_parameters','nLLVec','runlist_completed')
    end
end


%% GET ML PARAMETER ESTIMATES FOR EACH SUBJECT
clear all
filepath = 'fits/grace/';

subjVec = 1:9;
model = 'proportional';
TMScond = 1;
hemifield = 0;

switch model
    case 'proportional'
        nParams = 2;
    case 'min_error'
        nParams = 3;
end

nSubj = length(subjVec);

nll = nan(1,nSubj);
bfp = nan(nSubj,nParams);
for isubj = 1:nSubj;
    filename = [filepath 'fits_model_' num2str(model) '_subj' num2str(isubj) ...
        '_TMScond_' num2str(TMScond) 'hemi' num2str(hemifield) '.mat'];
    load(filename)
    
    blah = min(nLLVec);
    if length(blah) > 1; blah = blah(1); end
    idx = find(nLLVec == blah);
    nll(isubj) = blah;
    bfp(isubj,:) = ML_parameters(idx,:);
end

ML_parameters = bfp;
nLLVec = nll;
save([filepath 'fits_model_' num2str(model) '_TMScond_' num2str(TMScond) ...
    '_hemi' num2str(hemifield) '.mat'],'ML_parameters','nLLVec','TMScond','hemifield','subjVec','model')

%% PLOTS OF PARAMETER FITS!