

%% set up and save data in minimal useable format

clear all

% all_TMS_cond: 1: no tms, 2: ips2 3: l_spcs

subjnumVec = [1:9 11 12 16 17];
exclVec = [13 20 21 22];
nSubj = length(subjnumVec);

condition = 'noTMS';
hemifield = 2; % 1: left, 2: right

data = cell(1,nSubj);
nTrials = nan(nSubj,2);

for isubj = 1:nSubj
    subjnum = subjnumVec(isubj)
    
    load(sprintf('/data/TMS_Priority/TMS_Priority_behav_73019/subj%02d_%s_behav.mat',subjnum,condition))
    
    subjid = sprintf('S%02d',subjnum);
    
    % EXCLUDE TRIALS
    % based on preset exclusion criteria
    blah = cell2mat(cellfun(@(x) ~any(ismember(x,exclVec)),s_all.excl_trial,'UniformOutput',false));
    
    blah(s_all.i_sacc_rt<0.1 | s_all.i_sacc_rt>0.7) = 0; % reaction time exclusion
    blah(s_all.f_sacc_err>10) = 0;      % final saccade error
    blah(s_all.i_sacc_err>10) = 0;      % initial saccade error exclusion
    
    % additional subject/run/trial specific dropping (ask grace if need more clarification here)
    % r_num: run number. t_num: trial number
    switch condition
        case 'l_ips2'
            blah((s_all.r_num==2) & (subjnum==3) & ismember(s_all.t_num, [26 27 28 29 30 31  32  33  34 35 36])) = 0; %exclude run01 trials 26-36 due daq err (ips2 run02,sess1)
            blah((s_all.r_num==7) & (subjnum==6)) = 0;
        case 'l_spcs'
            blah((s_all.r_num==7) & (subjnum==2)) = 0;
            blah((s_all.r_num==1) & (subjnum==9)) = 0;  % exclude run01 due to coil slip (spcs sess 1)
            blah((s_all.r_num==1) & (subjnum==3) & ismember(s_all.t_num, [27 28 29 30 31  32  33  34  35 36])) = 0; % exclude run02 trials 27-36 due daq err (spcs run01,sess1)
            blah((s_all.r_num==1) & (subjnum==5)) = 0;  % exclude run01 due to coil "light" (spcs sess 1)
            blah((s_all.r_num==1) & (subjnum==6) & ismember(s_all.t_num, [33 34 35 36])) = 0;
    end
        
    % indices of high and low priority trials
    highpri_trials = (s_all.trialinfo(:,1) == 31) & blah;
    lowpri_trials = (s_all.trialinfo(:,1) == 32) & blah;
    
    % additional indexing for hemifield presentation
    if (hemifield)
        highpri_trials = highpri_trials & (s_all.trialinfo(:,2) == hemifield);
        lowpri_trials = lowpri_trials & (s_all.trialinfo(:,2) == hemifield);
    end
    
    % number of trials
    nTrials(isubj,1) = sum(highpri_trials);
    nTrials(isubj,2) = sum(lowpri_trials);
    
    data{isubj}.subjid = subjid;
    data{isubj}.use_trial = blah;
    data{isubj}.i_sacc_err = {s_all.i_sacc_err(highpri_trials) s_all.i_sacc_err(lowpri_trials)};
    data{isubj}.f_sacc_err = {s_all.f_sacc_err(highpri_trials) s_all.f_sacc_err(lowpri_trials)};
end

if (hemifield)
    save(sprintf('data/tms/data_allsubj_%s_hemifield%d.mat',condition,hemifield),'data','nTrials')
else
    save(sprintf('data/tms/data_allsubj_%s.mat',condition),'data','nTrials')  
end

%% GET DATA IN CORRECT FORMAT FOR FITTING

% % error
% % error_i = all_data.s_all.i_sacc_err;
% error_f = all_data.s_all.f_sacc_err;
% 
% % indicate which subject, condition, and hemifield to use
% subjnum = 2; % 1-9
% TMScond = 1; % 1 = no tms, 2 = ips2, 3 = spcs
% hemifield = 0; % 1: left hemi. 2: right hemi. 0: both hemi
% 
% % get data from a particular subject and condition (and hemifield?)
% idx = all_data.use_trial;
% idx = idx & (all_data.subj_all == subjnum);
% idx = idx & (all_data.TMS_cond_all == TMScond);
% if (hemifield); idx = idx & (all_data.s_all.trialinfo(:,2) == hemifield); end
% 
% data = cell(1,2);
% idxx = idx & (all_data.s_all.trialinfo(:,1) == 31); % high priority
% data{1} = error_f(idxx,:);
% idxx = idx & (all_data.s_all.trialinfo(:,1) == 32); % low priority
% data{2} = error_f(idxx,:);
% 
% % delete any nans
% data{1}(isnan(data{1})) = [];
% data{2}(isnan(data{2})) = [];
% 
% data

%% PLOT DATA

clear all

condition = 'noTMS';
errortype = 'f_sacc_err';
hemifield = 0;



% load data
if (hemifield)
    load(sprintf('data/tms/data_allsubj_%s_hemifield%d.mat',condition,hemifield))
else
    load(sprintf('data/tms/data_allsubj_%s.mat',condition))
end

nPriorities = 2;
nSubj = length(data);


% histograms per subjects
xlims = linspace(0,10,16);

colorMat = [1 0 0; 0 0 1];

error = cell(1,nPriorities);
error_M = nan(nSubj,nPriorities);
for isubj = 1:nSubj    
    
    for ipriority = 1:nPriorities
        
        % histogram of euclidean error
        datacounts = hist(data{isubj}.f_sacc_err{ipriority},xlims);
        error{ipriority}(isubj,:) = datacounts./sum(datacounts);
        error_M(isubj,ipriority) = mean(data{isubj}.f_sacc_err{ipriority});
    end
end

mean(error_M)

% =========== group plot =====================
meanerror = cellfun(@mean,error,'UniformOutput',false);
semerror = cellfun(@(x) std(x)./sqrt(size(x,1)),error,'UniformOutput',false);
figure; hold on
set(gcf,'Position',[28 504 500 350])

for ipriority = 1:nPriorities
    
    errorbar(xlims,meanerror{ipriority},semerror{ipriority},'Color',colorMat(ipriority,:),'LineStyle','none','LineWidth',1);
    
end
defaultplot
axis([0 10 0 0.5])
xlabel('error','FontSize',16)
set(gca,'YTick',[0:0.1:0.5],'FontSize',12);
ylabel('proportion','FontSize',16)


%% =====================================================================
%               MODELING
% ======================================================================

% - fit model
% - get bfp
% - plot fits

%% FIT MODEL

clear all

condition = 'l_spcs';
errortype = 'f_sacc_err';
model = 'flexible';
hemifield = 2;

subjnumVec = [1:9 11 12 16 17];
nSubj = length(subjnumVec);
runlist = 1:20;
runmax = 20;
exppriorityVec = [2/3 1/3];
fixparams = [];

% load data
if (hemifield)
    load(sprintf('data/tms/data_allsubj_%s_hemifield%d.mat',condition,hemifield))
else
    load(sprintf('data/tms/data_allsubj_%s.mat',condition))
end
    
for isubj = 1:nSubj
    
    % define current data
    datt = data{isubj}.(errortype);
    
    % file saving name
    subjnum = subjnumVec(isubj);
    if (hemifield)
        filename = sprintf('fits/tms/fits_model_%s_%s_hemifield%d_subjS%02d.mat',model,condition,hemifield,subjnum)
    else
        filename = sprintf('fits/tms/fits_model_%s_%s_subjS%02d.mat',model,condition,subjnum)
    end
    
    try load(filename); catch; ML_parameters = []; nLLVec = []; end
    try runlist_completed*2; catch; runlist_completed = []; end % seeing if runlist_completed exists yet
    
    
    for irun = 1:length(runlist);
        runlistt = runlist(irun);
        rng(runlistt)
        try
            [bfp,fval,rlc] = fit_parameters(model,datt,exppriorityVec,runlistt,runmax,fixparams);
            
            ML_parameters = [ML_parameters; bfp];
            nLLVec = [nLLVec fval];
            runlist_completed = [runlist_completed rlc];
            save(filename,'ML_parameters','nLLVec','runlist_completed')
        end
    end
end

%% GET ML PARAMETER ESTIMATES FOR EACH SUBJECT
clear all

subjnumVec = [1:9 11 12 16 17];
nSubj = length(subjnumVec);
model = 'flexible';
condition = 'l_ips2';
hemifield = 2;

switch model
    case 'proportional'
        nParams = 2;
    case {'min_error', 'flexible'}
        nParams = 3;
end

if (hemifield)
    filemod = sprintf('fits/tms/fits_model_%s_%s_hemifield%d',model,condition,hemifield);
else
    filemod = sprintf('fits/tms/fits_model_%s_%s',model,condition);
end
nll = nan(1,nSubj);
bfp = nan(nSubj,nParams);
for isubj = 1:nSubj;
    subjnum = subjnumVec(isubj);
    
    load(sprintf('%s_subjS%02d.mat',filemod,subjnum))
    
    blah = min(nLLVec);
    if length(blah) > 1; blah = blah(1); end
    idx = find(nLLVec == blah);
    nll(isubj) = blah;
    bfp(isubj,:) = ML_parameters(idx,:);
end

ML_parameters = bfp;
nLLVec = nll;

save(sprintf('%s.mat',filemod),'ML_parameters','nLLVec','condition','subjnumVec','model')


%% PLOTS OF PARAMETER FITS!

clear all; close all

% ======== LOAD DATA =======

model = 'flexible';
condition = 'l_spcs';
hemifield = 2; % 1: left hemi. 2: right hemi. 0: both hemi

loadpreddata = 0;
indvlplot = 1;
% nTrials = [100 50];
exppriorityVec = [2/3 1/3];

% load data and get in correct format
if (hemifield)
    load(sprintf('data/tms/data_allsubj_%s_hemifield%d.mat',condition,hemifield))
else
    load(sprintf('data/tms/data_allsubj_%s.mat',condition))
end
nSubj = length(data);

% error_f = all_data.s_all.f_sacc_err; % final saccade error

nPriorities = length(exppriorityVec);
nTrials = 1e3*ones(1,nPriorities); % how many trials to simulate per priority

% ====== get ML parameter estimate for isubj =======
if (hemifield)
    load(sprintf('fits/tms/fits_model_%s_%s_hemifield%d.mat',model,condition,hemifield))
else
    load(sprintf('fits/tms/fits_model_%s_%s.mat',model,condition))
end

% ====== SIMULATE DATA FOR PARTICIPANTS ====== 
filename = sprintf('fits/tms/modelpred_model_%s_%s.mat',model,condition);
% filename = [filepath 'modelpred_model_' model '_TMScond_' ...
%     num2str(TMScond) '_hemi' num2str(hemifield) '.mat'];
if (loadpreddata)
    load(filename,'preddata')
else
    preddata = simulate_data(model,1,ML_parameters,nTrials,exppriorityVec);
    save(filename,'preddata')
end

% histograms per subjects
xlims = linspace(0,10,16);

colorMat = [1 0 0; 0 0 1];

[error, simerror] = deal(cell(1,nPriorities));
if (indvlplot); figure; end;
for isubj = 1:nSubj
   
    subplot(5,3,isubj);
    
    for ipriority = 1:nPriorities
        
        % histogram of euclidean error
        datacounts = hist(data{isubj}.f_sacc_err{ipriority},xlims);
        simdatacounts = hist(preddata{isubj}{ipriority}(:,1),xlims);
        error{ipriority}(isubj,:) = datacounts./sum(datacounts);
        simerror{ipriority}(isubj,:) = simdatacounts./sum(simdatacounts);
        
        if (indvlplot)
            plot(xlims,error{ipriority}(isubj,:),'Color',colorMat(ipriority,:));
            hold on;
            plot(xlims,simerror{ipriority}(isubj,:),':','Color',colorMat(ipriority,:));
            defaultplot
            if mod(isubj,1) == 1; xlabel('euclidean error'); end
        end
        
    end
end

% =========== group plot =====================

figure;
colorMat = {'r','b','k'};

ha = tight_subplot(1,nPriorities,.03,[.26 .05],[.11 .05]);
set(gcf,'Position',[28 504 800 236])

meanerror = cellfun(@mean,error,'UniformOutput',false);
semerror = cellfun(@(x) std(x)./sqrt(size(x,1)),error,'UniformOutput',false);
meansimerror = cellfun(@mean,simerror,'UniformOutput',false);
semsimerror = cellfun(@(x) std(x)./sqrt(size(x,1)),simerror,'UniformOutput',false);

for ipriority = 1:nPriorities
    
    axes(ha(ipriority))
    
    fill([xlims fliplr(xlims)],[meansimerror{ipriority}-semsimerror{ipriority}...
        fliplr(meansimerror{ipriority}+semsimerror{ipriority})],colorMat{ipriority},'EdgeColor','none','FaceAlpha',0.4);
    hold on;
    errorbar(xlims,meanerror{ipriority},semerror{ipriority},'Color','k','LineStyle','none','LineWidth',1);
    defaultplot
    axis([0 10 0 0.5])
    
    xlabel('error','FontSize',16)
    set(ha(ipriority),'YTick',[0:0.1:0.5],'FontSize',12);
    if ipriority ~= 1
        set(ha(ipriority),'YTickLabel','');
    else
        ylabel('proportion','FontSize',16)
    end
end

