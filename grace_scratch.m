

%%

clear all

root = '/share/data/TMS_Priority';

% EXCLUDING 7,8,9,11,12  % FINAL SACCADE ANALYSIS 
%  subj = [1 2 3 4 5 6 16 17 20 21 22 24];
%  subjects = {'1','2','3','4','5','6','16','17','20','21','22','24'};

% FOR INITIAL SACCADE ANALYSIS % %%%%%%%%
% subj = [1 2 3 4 5 7 11 12 16 17 20 21 22 24];
% subjects = {'1','2','3','4','5','7','11','12','16','17','20','21','22','24'};

% ALL SUBJECTS (COMBINING INITIAL AND FINAL)
% allsubjidVec = [1:7 11 12 16 17 20 21 22 24]
subj = [1 2 3 4 5 6 7 11 12 16 17 20 21 22 24];
subjects = {'1','2','3','4','5','6','7','11','12','16','17','20','21','22','24'};

TMS_cond = {'noTMS','l_ips2','l_spcs'};  %which order should we load, plot conditions?


% TRIAL EXCLUSION CODES: 
WHICH_EXCL = [13 20];

% for now, let's use cat_struct to load/concatenate all data...
all_subj = nan(36*10*length(subj)*length(TMS_cond),1);
all_TMS_cond = nan(36*10*length(subj)*length(TMS_cond),1);

all_data = [];

startidx = 1;

for ss = 1:length(subj)
    
    for condidx = 1:length(TMS_cond)
        fn = sprintf('%s/TMS_Priority_behav_051721/subj%02.f_%s_behav.mat',root,subj(ss),TMS_cond{condidx});
        fprintf('Loading %s\n',fn);
        this_data = load(fn);
        
        this_subj = subj(ss);
        
        all_data = cat_struct(all_data,this_data);
        all_subj(startidx:(startidx-1+size(this_data.s_all.trialinfo,1))) = this_subj;
        all_TMS_cond(startidx:(startidx-1+size(this_data.s_all.trialinfo,1))) = condidx;
        
        startidx = startidx+size(this_data.s_all.trialinfo,1);
        
        clear this_subj this_data;
    end
end

% let's try this pattern for now
all_subj = all_subj(1:(startidx-1));
all_TMS_cond = all_TMS_cond(1:(startidx-1));
all_data.subj_all = all_subj;
all_data.TMS_cond_all = all_TMS_cond;

% determine which trials to include

all_data.use_trial = ~cellfun( @any, cellfun( @(a) ismember(a, WHICH_EXCL), all_data.s_all.excl_trial, 'UniformOutput',false));

% drop trials with very short (< 100 ms) or very long RT (> 1 s)
all_data.use_trial(all_data.s_all.i_sacc_rt<0.1 | all_data.s_all.i_sacc_rt>0.7) = 0;

all_data.use_trial(all_data.s_all.f_sacc_err>10) = 0;
all_data.use_trial(all_data.s_all.i_sacc_err>10) = 0;

all_data.use_trial(ismember(all_data.s_all.r_num, [7]) & all_TMS_cond == 2 & all_subj==6) = 0;
all_data.use_trial(all_data.s_all.r_num ==7 & all_TMS_cond == 3 & all_subj==2) = 0;
all_data.use_trial(all_data.s_all.r_num ==1 & all_TMS_cond == 3 & all_subj==9) = 0;%exclude run01 due to coil slip (spcs sess 1)
all_data.use_trial(all_data.s_all.r_num ==1 & all_TMS_cond == 3 & all_subj==3 & ismember(all_data.s_all.t_num, [27 28 29 30 31  32  33  34  35  36])) = 0; %exclude run02 trials 27-36 due daq err (spcs run01,sess1)
all_data.use_trial(all_data.s_all.r_num ==2 & all_TMS_cond == 2 & all_subj==3 & ismember(all_data.s_all.t_num, [26 27 28 29 30 31  32  33  34 35  36])) = 0; %exclude run01 trials 26-36 due daq err (ips2 run02,sess1)
all_data.use_trial(all_data.s_all.r_num ==1 & all_TMS_cond == 3 & all_subj==5) = 0; %exclude run01 due to coil "light" (spcs sess 1)
all_data.use_trial(all_data.s_all.r_num ==1 & all_TMS_cond == 3 & all_subj==6 & ismember(all_data.s_all.t_num, [33 34 35 36])) = 0;

% make them separate
allsubjidVec = [1:7 11 12 16 17 20 21 22 24];
subjid_max = max(allsubjidVec);
usetrial_cell = cell(1,subjid_max);
nSubjs = length(allsubjidVec);
condVec = {'noTMS','IPS2','sPCS'};
for isubj = 1:nSubjs
    subjid = allsubjidVec(isubj);
    idx_usetrial = (all_data.subj_all == subjid) & (all_data.use_trial); % useable trials for current subject
    
    for icond = 1:3 % for each tms condition
        
        cond = condVec{icond};
        idx_cond = (all_data.TMS_cond_all(idx_usetrial) == icond) & (all_data.s_all.trialinfo(idx_usetrial,2) == 2); % conditino & right visual hemi
        priVec = all_data.s_all.trialinfo(idx_cond,1);
        i_sacc_err = all_data.s_all.i_sacc_err(idx_cond);
        f_sacc_err = all_data.s_all.f_sacc_err(idx_cond);
        
        for ipri = 1:2 % for each priority level
            idx_pri = priVec == (30+ipri);
            
            blah = i_sacc_err(idx_pri);
            blah(isnan(blah)) = [];
            data.initial.(cond){ipri} = blah;
            
            blah = f_sacc_err(idx_pri);
            blah(isnan(blah)) = [];
            data.final.(cond){ipri} = blah;
        end
    end
    
    save(sprintf('data/tms/data_subjid%d.mat',subjid),'data')
end

%% things saved across stuff

clear all

condVec = {'noTMS','IPS2','sPCS'};
nConds = length(condVec);

% EXCLUDING 7,8,9,11,12  % FINAL SACCADE ANALYSIS 
subjidVec.final = [1 2 3 4 5 6 16 17 20 21 22 24];

% FOR INITIAL SACCADE ANALYSIS % %%%%%%%%
subjidVec.initial = [1 2 3 4 5 7 11 12 16 17 20 21 22 24];
nSubjs = structfun(@(x) length(x),subjidVec,'UniformOutput',false);

save('fittingsettings.mat')

% ALL SUBJECTS (COMBINING INITIAL AND FINAL)
% % allsubjidVec = [1:7 11 12 16 17 20 21 22 24]
% subj = [1 2 3 4 5 6 7 11 12 16 17 20 21 22 24];

%% set up and save data in minimal useable format

% clear all
% 
% % all_TMS_cond: 1: no tms, 2: ips2 3: l_spcs
% 
% subjnumVec = [1:9 11 12 16 17];
% exclVec = [13 20 21 22];
% nSubj = length(subjnumVec);
% 
% condition = 'noTMS';
% hemifield = 1; % 1: left, 2: right
% 
% data = cell(1,nSubj);
% nTrials = nan(nSubj,2);
% 
% for isubj = 1:nSubj
%     subjnum = subjnumVec(isubj)
%     
%     load(sprintf('/data/TMS_Priority/TMS_Priority_behav_73019/subj%02d_%s_behav.mat',subjnum,condition))
%     
%     subjid = sprintf('S%02d',subjnum);
%     
%     % EXCLUDE TRIALS
%     % based on preset exclusion criteria
%     blah = cell2mat(cellfun(@(x) ~any(ismember(x,exclVec)),s_all.excl_trial,'UniformOutput',false));
%     
%     blah(s_all.i_sacc_rt<0.1 | s_all.i_sacc_rt>0.7) = 0; % reaction time exclusion
%     blah(s_all.f_sacc_err>10) = 0;      % final saccade error
%     blah(s_all.i_sacc_err>10) = 0;      % initial saccade error exclusion
%     
%     % additional subject/run/trial specific dropping (ask grace if need more clarification here)
%     % r_num: run number. t_num: trial number
%     switch condition
%         case 'l_ips2'
%             blah((s_all.r_num==2) & (subjnum==3) & ismember(s_all.t_num, [26 27 28 29 30 31  32  33  34 35 36])) = 0; %exclude run01 trials 26-36 due daq err (ips2 run02,sess1)
%             blah((s_all.r_num==7) & (subjnum==6)) = 0;
%         case 'l_spcs'
%             blah((s_all.r_num==7) & (subjnum==2)) = 0;
%             blah((s_all.r_num==1) & (subjnum==9)) = 0;  % exclude run01 due to coil slip (spcs sess 1)
%             blah((s_all.r_num==1) & (subjnum==3) & ismember(s_all.t_num, [27 28 29 30 31  32  33  34  35 36])) = 0; % exclude run02 trials 27-36 due daq err (spcs run01,sess1)
%             blah((s_all.r_num==1) & (subjnum==5)) = 0;  % exclude run01 due to coil "light" (spcs sess 1)
%             blah((s_all.r_num==1) & (subjnum==6) & ismember(s_all.t_num, [33 34 35 36])) = 0;
%     end
%         
%     % indices of high and low priority trials
%     highpri_trials = (s_all.trialinfo(:,1) == 31) & blah;
%     lowpri_trials = (s_all.trialinfo(:,1) == 32) & blah;
%     
%     % additional indexing for hemifield presentation
%     if (hemifield)
%         highpri_trials = highpri_trials & (s_all.trialinfo(:,2) == hemifield);
%         lowpri_trials = lowpri_trials & (s_all.trialinfo(:,2) == hemifield);
%     end
%     
%     % number of trials
%     nTrials(isubj,1) = sum(highpri_trials);
%     nTrials(isubj,2) = sum(lowpri_trials);
%     
%     data{isubj}.subjid = subjid;
%     data{isubj}.use_trial = blah;
%     data{isubj}.i_sacc_err = {s_all.i_sacc_err(highpri_trials) s_all.i_sacc_err(lowpri_trials)};
%     data{isubj}.f_sacc_err = {s_all.f_sacc_err(highpri_trials) s_all.f_sacc_err(lowpri_trials)};
% end
% 
% if (hemifield)
%     save(sprintf('data/tms/data_allsubj_%s_hemifield%d.mat',condition,hemifield),'data','nTrials')
% else
%     save(sprintf('data/tms/data_allsubj_%s.mat',condition),'data','nTrials')  
% end

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
m_data = cellfun(@mean,error,'UniformOutput',false);
sem_data = cellfun(@(x) std(x)./sqrt(size(x,1)),error,'UniformOutput',false);
figure; hold on
set(gcf,'Position',[28 504 500 350])

for ipriority = 1:nPriorities
    
    errorbar(xlims,m_data{ipriority},sem_data{ipriority},'Color',colorMat(ipriority,:),'LineStyle','none','LineWidth',1);
    
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
load('fittingsettings.mat')
exptype = 'initial';
subjidVec = subjidVec.(exptype);
nSubjs = nSubjs.(exptype);

condition = 'noTMS';
% errortype = 'f_sacc_err';
model = 'flexible';
hemifield = 0;

runlist = 1:20;
runmax = 20;
exppriorityVec = [2/3 1/3];
fixparams = [];

% % load data
% if (hemifield)
%     load(sprintf('data/tms/data_allsubj_%s_hemifield%d.mat',condition,hemifield))
% else
%     load(sprintf('data/tms/data_allsubj_%s.mat',condition))
% end
    
for isubj = 1:nSubjs
    subjid = subjidVec(isubj);
    
    % load fitting data
    load(sprintf('data/tms/data_subjid%d.mat',subjid))
    data = data.(exptype).(condition);
    
    % file saving name
    if (hemifield)
        filename = sprintf('fits/tms/fits_model_%s_%s_hemifield%d_subj%d.mat',model,condition,hemifield,subjid)
    else
        filename = sprintf('fits/tms/fits_%s_model_%s_%s_subj%d.mat',exptype,model,condition,subjid);
    end
    
    try load(filename); catch; ML_parameters = []; nLLVec = []; end
    try runlist_completed*2; catch; runlist_completed = []; end % seeing if runlist_completed exists yet
    
    
    for irun = 1:length(runlist);
        runlistt = runlist(irun);
        rng(runlistt)
        try
            [bfp,fval,rlc] = fit_parameters(model,data,exppriorityVec,runlistt,runmax,fixparams);
            
            ML_parameters = [ML_parameters; bfp];
            nLLVec = [nLLVec fval];
            runlist_completed = [runlist_completed rlc];
            save(filename,'ML_parameters','nLLVec','runlist_completed')
        end
    end
end


%% fit stronger hypotheis: "flexible" model but w constraints
clear all

load('fittingsettings.mat')
exptype = 'initial';

runlist = 1:20;
runmax = 20;
exppriorityVec = [2/3 1/3];

condhypVec = [1 2 3];
fixparams = [];

for isubj = 1:nSubjs.(exptype)
    subjid = subjidVec.(exptype)(isubj);
    
    % load fitting data
    load(sprintf('data/tms/data_subjid%d.mat',subjid))
    data = data.(exptype);
    
    % file saving name
    filename = sprintf('fits/tms/fits_%s_model_strtonghypflex_condhyp%d%d%d_subj%d.mat',exptype,condhypVec(1),condhypVec(2),condhypVec(3),subjid);
    
    try load(filename); catch; ML_parameters = []; nLLVec = []; end
    try runlist_completed*2; catch; runlist_completed = []; end % seeing if runlist_completed exists yet
    
    
    for irun = 1:length(runlist);
        runlistt = runlist(irun);
        rng(runlistt)
        try
%             [bfp,fval,rlc] = fit_parameters(model,datt,exppriorityVec,runlistt,runmax,fixparams);
            [bfp,fval,rlc] = fit_tms_stronghyp(data,exppriorityVec,runlistt,runmax,fixparams,condhypVec);
            
            ML_parameters = [ML_parameters; bfp];
            nLLVec = [nLLVec fval];
            runlist_completed = [runlist_completed rlc];
            save(filename,'ML_parameters','nLLVec','runlist_completed')
        end
    end
end


%% GET ML PARAMETER ESTIMATES FOR EACH SUBJECT
clear all

load('fittingsettings.mat')
exptype = 'final';
nSubjs = nSubjs.(exptype);
subjidVec = subjidVec.(exptype);
% model = 'strtonghypflex';
model = 'flexible';
condition = 'noTMS';
condhypVec = [1 3 2];


switch model
    case 'proportional'
        nParams = 2;
    case {'min_error', 'flexible'}
        nParams = 3;
    case 'strtonghypflex'
        nParams = 5;
end

switch model
    case 'strtonghypflex'
        filemod = sprintf('fits/tms/fits_%s_model_%s_condhyp%d%d%d',exptype,model,condhypVec(1),condhypVec(2),condhypVec(3));
    otherwise
        filemod = sprintf('fits/tms/fits_%s_model_%s_%s',exptype,model,condition);
end
% if (hemifield)
%     filemod = sprintf('fits/tms/fits_%s_model_%s_%s_hemifield%d',exptype,model,condition,hemifield);
% else
%     filemod = sprintf('fits/tms/fits_%s_model_%s_%s',exptype,model,condition);
% end
nll = nan(1,nSubjs);
bfp = nan(nSubjs,nParams);
for isubj = 1:nSubjs;
    subjid = subjidVec(isubj);
    
    load(sprintf('%s_subj%d.mat',filemod,subjid))
    
    blah = min(nLLVec);
    if length(blah) > 1; blah = blah(1); end
    idx = find(nLLVec == blah);
    nll(isubj) = blah;
    bfp(isubj,:) = ML_parameters(idx,:);
end

ML_parameters = bfp;
nLLVec = nll;

save(sprintf('%s.mat',filemod),'ML_parameters','nLLVec','subjidVec','model')


%% PLOTS OF PARAMETER FITS: STRONG HYPOTHESIS

clear all

load('fittingsettings.mat')
exptype = 'initial';
nSubjs = nSubjs.(exptype);
subjidVec = subjidVec.(exptype);
model = 'strtonghypflex';
condhypVec = [1 3 2];

loadpreddata = 0;
indvlplot = 0;
exppriorityVec = [2/3 1/3];
nPriorities = length(exppriorityVec);
nTrials = [100 50];

load(sprintf('fits/tms/fits_%s_model_%s_condhyp%d%d%d.mat',exptype,model,condhypVec(1),condhypVec(2),condhypVec(3)));

xlims = linspace(0,10,8);
colorMat = {'g','b'};

[dMat.noTMS, dMat.IPS2, dMat.sPCS, mMat.noTMS, mMat.IPS2, mMat.sPCS,...
    datacounts.noTMS, datacounts.IPS2, datacounts.sPCS] = deal([]);
for isubj = 1:nSubjs
    subjid = subjidVec(isubj);
    
    % ======== LOAD DATA =======
    load(sprintf('data/tms/data_subjid%d.mat',subjid))
    data = data.(exptype);
    
    for icond = 1:nConds
        cond = condVec{icond};
        dMat.(cond) = [dMat.(cond); cellfun(@mean,data.(cond),'UniformOutput',true)];
        
        tempdatacounts = cellfun(@(x) hist(x,xlims),data.(cond),'UniformOutput',false);
        for ipriority = 1:nPriorities
            datacounts.(cond){ipriority}(isubj,:) = tempdatacounts{ipriority}./sum(tempdatacounts{ipriority});
        end
    end
    
    % ===== SIMUALTE DATA ======
    xbest = ML_parameters(isubj,:);
    
    for icond = 1:nConds
        cond = condVec{icond};
        condhyp = condhypVec(icond);
        
        switch condhyp
            case 1
                theta = xbest([1 2 3]);
            case 2
                theta = xbest([1 2 5]); % diff pVec
            case 3
                theta = xbest([4 2 3]); % diff Jbar
        end
        
        xx = simulate_data('flexible',1,theta,nTrials,exppriorityVec);
        mMat.(cond) = [mMat.(cond); cellfun(@mean,xx{1},'UniformOutput',true)];
        
        tempmodelcounts = cellfun(@(x) hist(x,xlims),xx{1},'UniformOutput',false);
        for ipriority = 1:nPriorities
            modelcounts.(cond){ipriority}(isubj,:) = tempmodelcounts{ipriority}./sum(tempmodelcounts{ipriority});
        end
    end
end

figure;
set(gcf,'Position',[267 496 1223 302])
for icond = 1:nConds
    cond = condVec{icond};
    
    subplot(1,3,icond); hold on;
    for ipriority = 1:nPriorities
        
        m_data = mean(datacounts.(cond){ipriority});
        sem_data = std(datacounts.(cond){ipriority})/sqrt(nSubjs);
        m_model = mean(modelcounts.(cond){ipriority});
        sem_model = std(modelcounts.(cond){ipriority})/sqrt(nSubjs);
        
        fill([xlims fliplr(xlims)],[m_model-sem_model fliplr(m_model+sem_model)],colorMat{ipriority},'EdgeColor','none','FaceAlpha',0.4);
        hold on;
        errorbar(xlims,m_data,sem_data,'Color','k','LineStyle','-','LineWidth',1);
    end
    defaultplot
%     switch exptype
%         case 'initial'
%             ylim([0 0.6])
%         case 'final'
            ylim([0 0.7])
%     end
    if (icond == 1);
        ylabel('proportion')
    end
    xlabel('error')
end

% figure; 
% for icond = 1:nConds
%     cond = condVec{icond};
%     
%     m_data = mean(dMat.(cond));
%     sem_data = std(dMat.(cond))/sqrt(nSubjs);
%     m_model = mean(mMat.(cond));
%     sem_model = std(mMat.(cond))/sqrt(nSubjs);
%     
%     subplot(1,3,icond)
%     errorbar(m_data,sem_data,'k-');  hold on;
%     for ipriority = 1:nPriorities
%         priority = exppriorityVec(ipriority);
%         fill([0.7 1.3 1.3 0.7]+ipriority-1,m_model(ipriority)+[-1 -1 1 1]*sem_model(ipriority),'r','EdgeColor','none','FaceAlpha',0.4)
%     end
%     switch exptype
%         case 'initial'
%             ylim([1.5 2.5])
%         case 'final'
%             ylim([1 1.8])
%     end
%     title(cond)
%     defaultplot
% end


%% PLOTS OF PARAMETER FITS! ONE-CONDITION MODELS

clear all; close all

load('fittingsettings.mat')
exptype = 'final';
nSubjs = nSubjs.(exptype);
subjidec = subjidVec.(exptype);

model = 'flexible';
condition = 'noTMS';
hemifield = 0; % 1: left hemi. 2: right hemi. 0: both hemi

loadpreddata = 0;
indvlplot = 0;
exppriorityVec = [2/3 1/3];
nPriorities = length(exppriorityVec);
nTrials = 1e3*ones(1,nPriorities); % how many trials to simulate per priority


% plotssing stuff
xlims = linspace(0,10,16);
colorMat = [1 0 0; 0 0 1];

% load ML parameters
load(sprintf('fits/tms/fits_%s_model_%s_%s.mat',exptype,model,condition));

[dMat, mMat] = deal(cell(1,nPriorities));
for isubj = 1:nSubjs
    subjid = subjidVec(isubj);
    
    % ======== LOAD DATA =======
    load(sprintf('data/tms/data_subjid%d.mat',subjid))
    data = data.(exptype).(condition);
    
    
    
    % ======= simulate data ======
    theta = ML_parameters(isubj,:);
    preddata = simulate_data(model,1,theta,nTrials,exppriorityVec);
    
    
    for ipriority = 1:nPriorities
        % histogram of euclidean error
        datacounts = hist(data{ipriority},xlims);
        simdatacounts = hist(preddata{1}{ipriority}(:,1),xlims);
        dMat{ipriority}(isubj,:) = datacounts./sum(datacounts);
        mMat{ipriority}(isubj,:) = simdatacounts./sum(simdatacounts);
        
    end
end


% ======== LOAD DATA =======



% % load data and get in correct format
% if (hemifield)
%     load(sprintf('data/tms/zz_old/data_allsubj_%s_hemifield%d.mat',condition,hemifield))
% else
%     load(sprintf('data/tms/zz_old/data_allsubj_%s.mat',condition))
% end
% nSubj = length(data);

% nPriorities = length(exppriorityVec);
% nTrials = 1e3*ones(1,nPriorities); % how many trials to simulate per priority
% 
% % ====== get ML parameter estimate for isubj =======
% if (hemifield)
%     load(sprintf('fits/tms/zz_old2/fits_model_%s_%s_hemifield%d.mat',model,condition,hemifield))
% else
%     load(sprintf('fits/tms/zz_old2/fits_model_%s_%s.mat',model,condition))
% end

% % ====== SIMULATE DATA FOR PARTICIPANTS ====== 
% filename = sprintf('fits/tms/zz_old2/modelpred_model_%s_%s.mat',model,condition);
% % filename = [filepath 'modelpred_model_' model '_TMScond_' ...
% %     num2str(TMScond) '_hemi' num2str(hemifield) '.mat'];
% if (loadpreddata)
%     load(filename,'preddata')
% else
%     preddata = simulate_data(model,1,ML_parameters,nTrials,exppriorityVec);
%     save(filename,'preddata')
% end

% % histograms per subjects
% xlims = linspace(0,10,16);
% 
% colorMat = [1 0 0; 0 0 1];

% if (indvlplot); figure; end;
% for isubj = 1:nSubj
%    
%     subplot(5,3,isubj);
%     
%     for ipriority = 1:nPriorities
%         
%         % histogram of euclidean error
%         datacounts = hist(data{isubj}.f_sacc_err{ipriority},xlims);
%         simdatacounts = hist(preddata{isubj}{ipriority}(:,1),xlims);
%         error{ipriority}(isubj,:) = datacounts./sum(datacounts);
%         simerror{ipriority}(isubj,:) = simdatacounts./sum(simdatacounts);
%         
%         if (indvlplot)
%             plot(xlims,error{ipriority}(isubj,:),'Color',colorMat(ipriority,:));
%             hold on;
%             plot(xlims,simerror{ipriority}(isubj,:),':','Color',colorMat(ipriority,:));
%             defaultplot
%             if mod(isubj,1) == 1; xlabel('euclidean error'); end
%         end
%         
%     end
% end

% =========== group plot =====================

figure;
colorMat = {'r','b','k'};

ha = tight_subplot(1,nPriorities,.03,[.26 .05],[.11 .05]);
set(gcf,'Position',[28 504 800 236])

m_data = cellfun(@mean,dMat,'UniformOutput',false);
sem_data = cellfun(@(x) std(x)./sqrt(size(x,1)),dMat,'UniformOutput',false);
m_model = cellfun(@mean,mMat,'UniformOutput',false);
sem_model = cellfun(@(x) std(x)./sqrt(size(x,1)),mMat,'UniformOutput',false);

for ipriority = 1:nPriorities
    
    axes(ha(ipriority))
    
    fill([xlims fliplr(xlims)],[m_model{ipriority}-sem_model{ipriority}...
        fliplr(m_model{ipriority}+sem_model{ipriority})],colorMat{ipriority},'EdgeColor','none','FaceAlpha',0.4);
    hold on;
    errorbar(xlims,m_data{ipriority},sem_data{ipriority},'Color','k','LineStyle','none','LineWidth',1);
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

%% MODEL COMPARISON

clear all
exptype = 'final';

modelVec = {'strtonghypflex_condhyp123','strtonghypflex_condhyp132'};
nModels = length(modelVec);

for imodel = 1:nModels
    model = modelVec{imodel};
    load(sprintf('fits/tms/fits_%s_model_%s.mat',exptype,model));
    
    LLMat(imodel,:) = nLLVec;
end

LLMat = bsxfun(@minus,LLMat,LLMat(imodel,:));

figure;
plot(LLMat,'o')
defaultplot
xlim([0.7 2.3])
