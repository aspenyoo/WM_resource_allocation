%% CLUSTER RELATED

clear all
nSubj = 11;
expnumber = 1;
model = 3;

filepath = ['fits/exp' num2str(expnumber) '/'];

for isubj = 1:nSubj
    isubj;
    load([filepath 'fits_model' num2str(model) '_subj' num2str(isubj) '.mat'])
    
    blah = sort(unique(runlist_completed));
%     length(blah)
%     pause
    bleh(isubj) = length(blah);
end
bleh

%% % % % % % % % % % % % % % % % % % % % % 
%       DATA RELATED
% % % % % % % % % % % % % % % % % % % % %

%% looking at nTrials

for isubj = 1:11
    for ipriority = 1:3
        blah(isubj,ipriority) = size(data{isubj}{ipriority},1);
    end
end
blah
bleh = bsxfun(@rdivide,blah(:,1), blah)
mean(bleh)

%% get data in clean format for data fitting

% load in data file
% load('group_data.mat')

% saving subjid, priority, radius, response_x, response_y, target_theta
usefuldata = group_data(:,[1 2 11 12 13 14]);

% deleting nans
idx = logical(sum(isnan(usefuldata),2)); % idxs of nans
usefuldata(idx,:) = [];

% rename useful columns of group_data.mat
data_subj = usefuldata(:,1);
data_priority = usefuldata(:,2);
% data_radius = usefuldata(:,3);
response_x = usefuldata(:,3);
response_y = usefuldata(:,4);
target_x = usefuldata(:,5);
target_y = usefuldata(:,6);
% response_x = usefuldata(:,4);
% response_y = usefuldata(:,5);
% target_x = 10*cosd(usefuldata(:,6)); % x location
% target_y = 10*sind(usefuldata(:,6)); % y location

% errors
error_x = response_x - target_x;
error_y = response_y - target_y;
error_distance = sqrt(error_x.^2 + error_y.^2);

% deleting trials with technical malfunctions
idx = (error_distance >= 10);
% idx = (error_distance >= 10) | (data_radius >= 10);
error_distance(idx) = [];
data_priority(idx) = [];
data_subj(idx) = [];
% data_radius(idx) = [];
error_x(idx) = [];
error_y(idx) = [];

% priority and subject numbers
priorityVec = [0.6 0.3 0.1]; % sort(unique(data_priority),'descend'); % priority condition
nPriorities = length(priorityVec);
titleVec = {'high','med','low'};
subjVec = unique(data_subj);
nSubj = length(subjVec);

data = cell(1,nSubj);
for isubj = 1:nSubj
    subjnum = subjVec(isubj);
    
    data{isubj} = cell(1,nPriorities);
    for ipriority = 1:nPriorities
        priority = priorityVec(ipriority);
        idx = (data_subj == subjnum) & (data_priority == priority);
        nTrials = sum(idx);
        
        %         data{isubj}{ipriority} = nan(nTrials,2);
        data{isubj}{ipriority}(:,1) = error_distance(idx);
        %         data{isubj}{ipriority}(:,2) = data_radius(idx);
    end
end

save('cleandata_nodisc.mat','data')

%% check jitter

[t, r] = cart2pol(group_data(:,8),group_data(:,9));
plot(t,'o')


%% look at guess distributions. permutation test

clear all
isubj = 4;
expnumber = 2;
priorityVec = [0.1 0.3 0.6];
nPriorities = length(priorityVec);
nbins = 15;
xmax = 30;

load(['exp' num2str(expnumber) '_groupdata.mat'])


figure;
for ipriority = 1:nPriorities
    priority = priorityVec(ipriority);
    
    idx = (group_data(:,2) == priority) & (group_data(:,1) == isubj);
    
    switch expnumber
        case 1
            x_resp = group_data(idx,11);
            y_resp = group_data(idx,12);
            x_target = group_data(idx,13);
            y_target = group_data(idx,14);
        case 2
            x_resp = group_data(idx,6);
            y_resp = group_data(idx,7);
            [x_target,y_target] = pol2cart(group_data(idx,13)/180*pi,10);
%             x_target = group_data(idx,8);
%             y_target = group_data(idx,9);
    end
    
    % ========== actual data =========
    % calculate error of actual data
    error_x = x_resp - x_target;
    error_y = y_resp - y_target;
    error_distance = sqrt(error_x.^2 + error_y.^2);
    
    % plot distributions of errors
    subplot(2,3,ipriority)
    histogram(error_distance,linspace(0,xmax,nbins));
    xlim([0 xmax])
    
    % =========== permuted data ========
    % permute data
    idxx = randperm(sum(idx));
    x_respperm = x_resp(idxx);
    y_respperm = y_resp(idxx);
    
    % calculate error
    error_x = x_respperm - x_target;
    error_y = y_respperm - y_target;
    error_distance = sqrt(error_x.^2 + error_y.^2);
    
    % plot distribution of errors
    subplot(2,3,3+ipriority)
    histogram(error_distance,linspace(0,xmax,nbins));
    xlim([0 xmax])
    
end

%% bootstrap guess distribution for each priority

nBoots = 10;

counts = cell(1,3);
for ipriority = 1:nPriorities
    priority = priorityVec(ipriority)
    counts{ipriority} = nan(nBoots,nbins-1);
    
    for iboot = 1:nBoots
        
        % permute data
        idxx = randperm(sum(idx));
        x_respperm = x_resp(idxx);
        y_respperm = y_resp(idxx);
        
        % calculate error
        error_x = x_respperm - x_target;
        error_y = y_respperm - y_target;
        error_distance = sqrt(error_x.^2 + error_y.^2);
        
        % get histogram counts
        counts{ipriority}(iboot,:)  = histcounts(error_distance,linspace(0,xmax,nbins));
        
    end
    
    % plot M sem of counts
    M = mean(counts{ipriority});
    SEM = std(counts{ipriority})/sqrt(nBoots);
    centers = linspace(0,xmax,nbins);
    centers = centers(1:end-1) + diff(centers(1:2))/2;
    subplot(1,3,ipriority)
    plot_summaryfit(centers,[],[],M,SEM,aspencolors('booger'));

end

%% nTrials for each subject
clear all
expnumber = 1;

filename = ['exp' num2str(expnumber) '_cleandata.mat'];
load(filename)

nSubj = length(data);
nTrials = nan(nSubj,3);
for isubj = 1:nSubj
    for ipriority = 1:3
        nTrials(isubj,ipriority) = size(data{isubj}{ipriority},1);
    end
end

save(filename,'data','nTrials')
% [AIC,BIC,AICc] = modcomp(nLL,K,n);


%% make excel files into one big mat file
clear all
expnumber = 2;

% get all file names in this folder
filepath = ['exp' num2str(expnumber) '_rawdata/'];
blah = dir([filepath '*.csv']);
nFiles = length(blah);

concatMat = [];
for ifile = 1:nFiles
    filename = blah(ifile).name
    
    T = readtable([filepath filename],'ReadVariableNames',true);
    mat = table2array(T(:,1:end-1));
    
%     mat = dlmread(filename,',',1,0);
%     size(concatMat)
%     size(mat)
    concatMat = [concatMat; mat];
end

data = concatMat;
names = T.Properties.VariableNames(1:end-1);
% names = {'subject'	'run'	'trial'	'response'	'delay'	'RT'	'distance'	'target_value'	'target_loc_quadrant'	'target_loc_angle'	'distractor1_value'	'distractor1_loc_quadrant'	'distractor1_loc_angle'	'distractor2_value'	'distractor2_loc_quadrant'	'distractor2_loc_angle'	'target_value'	'distractor3_loc_quadrant'	'distractor3_loc_angle'	'order_0'	'order_0.1'	'order_0.3'	'order_0.6' 'timestamp'};
save(['exp' num2str(expnumber) '_rawdata.mat'],'data','names');

% names = {'subject'	'run'	'trial'	'accuracy'	'delay'	'disc_size'	'reward_at_stake'	'target_value'	'error'	'memorandum_x'	'memorandum_y'	'actual_location_x'	'actual_location_y'	'target_loc_quadrant'	'target_loc_angle'	'distractor1_value'	'distractor1_loc_quadrant'	'distractor1_loc_angle'	'distractor2_value'	'distractor2_loc_quadrant'	'distractor2_loc_angle'	'distractor3_value'	'distractor3_loc_quadrant'	'distractor3_loc_angle'};

%% check bias as a function of distractor priority

clear all

isubj = 9;
target = 0.3;
distractor = 0.6;

load('exp1_rawdata.mat')

% keep only isubj's data
idx = data(:,1) == isubj; 
data = data(idx,:); 

% % delete those with RT higher than some cutoff
% cutoff = 2500;
% idx = data(:,5) > cutoff; 
% data(idx) = [];

% get target loc (just angle)
idx = data(:,8) == target;
data = data(idx,:);
targetloc = data(:,10);

% get disctractor loc (just angle)
idx = data == distractor;
idx = [idx(:,end-1:end) idx(:,1:end-2)]';
dataT = data';
distractorloc = dataT(idx);

% response (angle)

%% 

load('exp1_zuzprocesseddata.mat')

idx = (group_data(:,1) == isubj) & (group_data(:,2) == target);

%% % % % % % % % % % % % % % % % % % % % % % % % % %
%       MODEL RELATED
% % % % % % % % % % % % % % % % % % % % % % % % % %

%% understanding Jbar, tau, and sigma

isubj = 1;
ipriority = 3;
Jbar = JbarMat(isubj,ipriority);
tau = ML_parameters(isubj,2);
beta = ML_parameters(isubj,end);

% Jbar = 1;
% tau = .1;
% beta = 0.01;

[JVec] = loadvar('JVec',{Jbar,tau});
% JVec = linspace(1,2,100);
Jpdf = gampdf(JVec,Jbar/tau,tau);
Jpdf = Jpdf./qtrapz(Jpdf); % normalize


% figure;
subplot(3,1,1) % JVec
plot(JVec,Jpdf,'k-'); defaultplot
subplot(3,1,2) % relative to mean, Jbar
plot(JVec./Jbar,Jpdf,'k-'); defaultplot
subplot(3,1,3) % sigma
plot(1./sqrt(JVec),Jpdf,'k-'); defaultplot


%% look at typical max for J/Jbar gamma distribution

expnumber = 2;
imodel = 1;

load(['fits/exp' num2str(expnumber) '/fits_model' num2str(imodel) '.mat'])

switch expnumber
    case 1
        nSubj = 14;
    case 2
        nSubj = 11;
end

JVecMax = nan(nSubj,3);
for isubj = 1:nSubj
    isubj
    
    Theta = ML_parameters(isubj,:);
    
    switch imodel
        case 1 % optimal
            % calculate the proportions that maximize expected utility
            pVec = calc_optimal_pVec(Theta);
        case 2 % not optimal
            if sum(Theta(end-1:end))>1 % reflect over pHigh + pMed = 1 line
                pVec = [1-Theta(end) 1-Theta(end-1)];
                pVec = [pVec 1-sum(pVec)];
            else
                pVec = [Theta(end-1:end) 1-sum(Theta(end-1:end))];
            end
        case 3 % fixed
            pVec = [0.6 0.3 0.1];
    end
    
    tau = Theta(2);
    JVec = linspace(0,10,100);
    for iJbar = 1:3
        Jbar = Theta(1)*pVec(iJbar);
        %         [JVec] = loadvar({'JVec',Jbar,tau});
        
        %         Jpdf = gampdf(JVec,Jbar/tau,tau);
        Jpdf = gampdf(JVec*Jbar,Jbar/tau,tau);
        Jpdf = Jpdf./qtrapz(Jpdf); % normalize
        
        %         JVecMax(isubj,iJbar) = JVec(end)/Jbar;
        
        plot(JVec,Jpdf,'k-'); defaultplot;
        pause;
        %         plot(JVec./Jbar,Jpdf,'k-'); defaultplot;
        %         pause;
    end
end

% JVecMax

%% optimal pVec as a function of Jbar_total
clear

tau = 1;
alpha = 1;
beta = 1;
N = 20;
JbartotalVec = linspace(1e-3,10,N);

pVec = nan(N,3);
for iJbartotal = 1:N
    iJbartotal
    
    Jbartotal = JbartotalVec(iJbartotal);
    pVec(iJbartotal,:) = calc_optimal_pVec([Jbartotal tau alpha beta]);
end

plot(bsxfun(@times,pVec,JbartotalVec'),'k-')

%%
clear all

Jbar = 2;
tau = 1;
alpha = 1;
beta = .1;

J = gamrnd(Jbar,Jbar/tau)
rVec = loadvar('rVec');

EU = calc_EU(rVec,J,alpha);
pdf_r = calc_pdf_r(beta,J,alpha);

subplot(2,1,1)
plot(rVec,EU)

subplot(2,1,2)
plot(rVec,pdf_r)


%% looking at expected number of points
clear all

Jbar = 2;
tau = 1;
alpha = 1;
beta = 0.5;

[JVec, rVec] = loadvar('JVec',{Jbar,tau},'rVec');
Jpdf = gampdf(JVec,Jbar/tau,tau);
Jpdf = Jpdf./qtrapz(Jpdf); % normalize
pdf_r = calc_pdf_r(beta, JVec, alpha);

EU = calc_EU(rVec,JVec,alpha);


%% % % % % % % % % % % % % % % % % % % % % % %
%       MODEL FIT RELATED
% % % % % % % % % % % % % % % % % % % % % % %

% - get ML parameter estimates
% - recalculate nLLs
% - calculate pVecs for different models
% - ternary plots


%% see what runlist idxs you need to do still
clear all; clc

expnumber = 1;
imodel = 4;
nSubj = 14;

filepath = ['fits/exp' num2str(expnumber) '/'];

for isubj = 1:nSubj
    load([filepath 'fits_model' num2str(imodel) '_subj' num2str(isubj) '.mat'])

    blah = 1:50;
    blah(runlist_completed) = [];
    fprintf([num2str(isubj) ' ' num2str(blah) '\r'])
end

%% check fits of indvl subjects

clear all

expnumber = 2;
imodel = 1;
isubj = 3;

load(sprintf('data/exp%d_cleandata.mat',expnumber))


load(sprintf('fits/exp%d/fits_model%d_subj%d.mat',expnumber,imodel,isubj))

nLLcol = size(ML_parameters,2)+1;
blah = sortrows([ML_parameters nLLVec'],nLLcol);
blah(:,1:2) = log(blah(:,1:2));
blah = [blah nan(size(blah,1),1)];

for iblah = 1:10; 
    blah(iblah,end) = calc_nLL(1,blah(iblah,1:4),data{isubj}); 
    blah(iblah,5:end)
end

blah

%%   GET ML PARAMETER ESTIMATES

clear all

expnumber = 2;
subjVec = 1:10;
imodel = 4;


testmodel = 4;
fakedata = 1;
isriskfixed = 0;
nSubj = length(subjVec);

if (isriskfixed)
    filepath = ['fits/exp' num2str(expnumber) '_fixedrisk/'];
else
    filepath = ['fits/exp' num2str(expnumber) '/'];
end

if (fakedata)
    pretxt = 'modelrecov';
else
    pretxt = 'fits';
end

switch testmodel
    case {1,3}
        nParams = 4;
    case 2
        nParams = 6;
    case 4
        nParams = 5;
end
if (expnumber == 1); nParams = nParams - 2; end

bfp = nan(nSubj,nParams);
nLL = nan(1,nSubj);
for subjnum = 1:nSubj
    isubj = subjVec(subjnum);
    
    if strcmp(pretxt,'modelrecov')
        load([filepath pretxt '_truemodel' num2str(imodel) '_testmodel' num2str(testmodel) '_subj' num2str(isubj) '.mat'],'ML_parameters','nLLVec')
    else
        load([filepath pretxt '_model' num2str(imodel) '_subj' num2str(isubj) '.mat'],'ML_parameters','nLLVec')
    end
    %     load([filepath pretxt '_model' num2str(imodel) '_subj' num2str(isubj) '.mat'])
    blah = ML_parameters(nLLVec == min(nLLVec),:);
    bfp(subjnum,:) = blah(1,:);
    nLL(subjnum) = min(nLLVec);
end

ML_parameters = bfp;
nLLVec = nLL;
if strcmp(pretxt,'modelrecov')
    save([filepath pretxt '_truemodel' num2str(imodel) '_testmodel' num2str(testmodel) '.mat'],'ML_parameters','nLLVec')
else
    save([filepath pretxt '_model' num2str(imodel) '.mat'],'ML_parameters','nLLVec')
end

%% fit exp1 

clear all
testmodel = 2;
truemodel = testmodel;
runmax = 50;
expnumber = 1;
switch expnumber
    case 1
        nSubj = 14;
    case 2
        nSubj = 11;
end


for isubj = 1:nSubj
    isubj
    
%     load(['fits/exp' num2str(expnumber) '/fits_model' num2str(testmodel) '_subj' num2str(isubj) '.mat'])
    runlist = 1:runmax;
%     runlist(unique(runlist_completed)) = [];
    
    fit_parameters(testmodel,isubj,runlist,runmax,truemodel,expnumber)
end


%% recalculate nLL for each subject

clear all
expnumber = 1;
imodel = 4;
nSubj = 14;

load(['exp' num2str(expnumber) '_cleandata.mat'])
filepath = ['fits/exp' num2str(expnumber) '/'];

for isubj = 1:nSubj
    isubj
    
    load([filepath 'fits_model' num2str(imodel) '_subj' num2str(isubj) '.mat'])
    ML_parameters(:,1:2) = log(ML_parameters(:,1:2));
    
    nStarts = length(nLLVec);
    newnLL = nan(1,nStarts);
    for inll = 1:nStarts
        fprintf([num2str(inll) ', '])
        newnLL(inll) = calc_nLL(imodel,ML_parameters(inll,:),data{isubj});
    end 
    
    [nLLVec; newnLL]
%     nLLVec = newnLL;
%     ML_parameters(:,1:2) = exp(ML_parameters(:,1:2));
%     save([filepath 'fits_model' num2str(imodel) '_subj' num2str(isubj) '.mat'],...
%         'runlist_completed','nLLVec','ML_parameters')
end

%% recalculate best fit nLL

clear all
expnumber = 1;
modelVec = [4];

load(['exp' num2str(expnumber) '_cleandata.mat'])
filepath = ['fits/exp' num2str(expnumber) '/'];

for imodel = modelVec
    load([filepath 'fits_model' num2str(imodel) '.mat'])
    ML_parameters(:,1:2) = log(ML_parameters(:,1:2));
  
    nSubj = length(nLLVec);
    newnLL = nan(1,nSubj);
    for isubj = 1:nSubj
        fprintf([num2str(isubj) ', '])
        newnLL(isubj) = calc_nLL(imodel,ML_parameters(isubj,:),data{isubj});
    end 

    [newnLL; nLLVec]
%     nLLVec = newnLL;
%     ML_parameters(:,1:2) = exp(ML_parameters(:,1:2));
%     save([filepath 'fits_model' num2str(imodel) '_subj' num2str(isubj) '.mat'],...
%         'nLLVec','ML_parameters')
end

%% check nLL

% clear all
expnumber = 2;
imodel = 1;
subjVec = 1:11;

nSubj = length(subjVec);
load(['exp' num2str(expnumber) '_cleandata.mat'])

filepath = ['fits/exp' num2str(expnumber) '/'];
load([filepath 'fits_model' num2str(imodel) '.mat'])
ML_parameters(:,1:2) = log(ML_parameters(:,1:2));

switch expnumber
    case 1
        nLL4 = nan(1,14);
    case 2
        nLL4 = nan(1,11);
end
for isubj = 1:nSubj
    subjnum = subjVec(isubj);
    subjnum
    
    nLL4(subjnum) = calc_nLL(imodel,ML_parameters(subjnum,:),data{subjnum});
end

[nLLVec; nLL4]

%% optimal pVec for model 1

clear all
imodel = 1;
expnumber = 2;

load(['fits/exp' num2str(expnumber) '/fits_model' num2str(imodel) '.mat'])
nSubj = 11;

pMat = nan(nSubj,3);
for isubj = [1 6 10 11]
    isubj
    pMat(isubj,:) = calc_optimal_pVec(ML_parameters(isubj,:));
end

save(['fits/exp' num2str(expnumber) '/fits_model' num2str(imodel) '.mat'],'ML_parameters','nLLVec')

figure;
% plot(pMat)
plot(bsxfun(@times,pMat,ML_parameters(:,1)))
defaultplot
xlabel('subject number')
ylabel('Jbar_{[condition]}')

%% pVec for model2

clear all
imodel = 2;
expnumber = 2;
fixedrisk = '';
filepath = ['fits/exp' num2str(expnumber)  fixedrisk '/'];
load([filepath 'fits_model' num2str(imodel) '.mat'])
nSubj = size(ML_parameters,1);

pMat = ML_parameters(:,end-1:end);
pMat(:,3) = 1-sum(pMat,2);

% plot(pMat)
figure
plot(bsxfun(@times,pMat,ML_parameters(:,1)))
hold on;

%% pVec for model4

clear all
imodel = 4;
expnumber = 1;
filepath = ['fits/exp' num2str(expnumber) '/'];
load([filepath 'fits_model' num2str(imodel) '.mat'])
nSubj = size(ML_parameters,1);

pMat = nan(nSubj,3);
for isubj = 1:nSubj;
    pMat(isubj,:) = calc_pVec_minerror(ML_parameters(isubj,:));
end

% plot(pMat)
figure
plot(bsxfun(@times,pMat,ML_parameters(:,1)))
hold on;


%% compare flexible model LL with parameters from other model 

load(['exp' num2str(expnumber) '_cleandata.mat'])
[logflag] = loadconstraints(imodel,expnumber);

ML_parameters(:,logflag) = log(ML_parameters(:,logflag));

for isubj = 1:14;
    isubj
    subjdata = data{isubj};
    theta = [ML_parameters(isubj,:) pMat(isubj,1:2)];
    
    flexnLLVec(2,isubj) = calc_nLL(2,theta,subjdata);
    realnLLVec(2,isubj) = calc_nLL(imodel,ML_parameters(isubj,:),subjdata);
end

%% plot triangle (ternary) plot with lines indicating 0.6 0.3 0.1

figure;
% axis
[h,hg,htick]=terplot;
c1 = [1 0.5 1/3];
c2 = [0 0.5 1/3];
c3 = [0 0 1/3];
hlabels=terlabel('high','medium','low');

% colors
orange = aspencolors('orange');
aqua = aspencolors('bluegreen');
indigo = aspencolors('indigo');

%
set(h,'LineWidth',1)

%--  Change the color of the grid lines
set(hg(:,1),'color',orange)
set(hg(:,2),'color',aqua)
set(hg(:,3),'color',indigo)

% make 0.6 0.3 0.1 solid
set(hg(3,1),'LineStyle','--')
set(hg(6,2),'LineStyle','--')
set(hg(1,3),'LineStyle','--')

%--  Modify the labels
set(hlabels,'fontsize',12)
set(hlabels(1),'color',orange)
set(hlabels(2),'color',aqua)
set(hlabels(3),'color',indigo)

%--  Modify the tick labels
set(htick(:,1),'color',orange,'linewidth',3)
set(htick(:,2),'color',aqua,'linewidth',3)
set(htick(:,3),'color',indigo,'linewidth',3)

% plot data
hter=ternaryc(pMat(:,1),pMat(:,2),pMat(:,3));
set(hter,'marker','o','markerfacecolor','k','markersize',4','markeredgecolor','k')

%% plot triangle (ternary) plot with monotonic area shaded

figure;
% axis
[h,hg,htick]=terplot;
c1 = [1 0.5 1/3];
c2 = [0 0.5 1/3];
c3 = [0 0 1/3];

x=0.5-c1*cos(pi/3)+c2/2;
y=0.866-c1*sin(pi/3)-c2*cot(pi/6)/2;
patch('Faces',[1 2 3],'Vertices',[x' y'],'FaceColor',aspencolors('seacolored'),'FaceAlpha',0.3,'EdgeColor','none');

% plot data
hter=ternaryc(pMat(:,1),pMat(:,2),pMat(:,3));
set(hter,'marker','o','markerfacecolor','k','markersize',4','markeredgecolor','k')
hold on;
hmod3 = ternaryc(0.6,0.3,0.1);
set(hmod3,'marker','o','markerfacecolor','r','markersize',4','markeredgecolor','r')
hlabels=terlabel('high','medium','low');

%% plot triangle (ternary) plot with three criteria triangles

figure;
% axis
[h,hg,htick]=terplot;

% high > med
c1 = [1 0.5 0];
c2 = [0 0.5 0];
c3 = [0 0 1];
x=0.5-c1*cos(pi/3)+c2/2;
y=0.866-c1*sin(pi/3)-c2*cot(pi/6)/2;
patch('Faces',[1 2 3],'Vertices',[x' y'],'FaceColor',[1 0 0],'FaceAlpha',0.3,'EdgeColor','none');

% med > low
c1 = [1 0 0];
c2 = [0 1 0.5];
c3 = [0 0 0.5];
x=0.5-c1*cos(pi/3)+c2/2;
y=0.866-c1*sin(pi/3)-c2*cot(pi/6)/2;
patch('Faces',[1 2 3],'Vertices',[x' y'],'FaceColor',[0 1 0],'FaceAlpha',0.3,'EdgeColor','none');

% high > low
c1 = [1 0 0.5];
c2 = [0 1 0];
c3 = [0 0 0.5];
x=0.5-c1*cos(pi/3)+c2/2;
y=0.866-c1*sin(pi/3)-c2*cot(pi/6)/2;
patch('Faces',[1 2 3],'Vertices',[x' y'],'FaceColor',[0 0 1],'FaceAlpha',0.3,'EdgeColor','none');

% plot data
hter=ternaryc(pMat(:,1),pMat(:,2),pMat(:,3));
set(hter,'marker','o','markerfacecolor','k','markersize',4','markeredgecolor','k')
hold on;
hmod3 = ternaryc(0.6,0.3,0.1);
set(hmod3,'marker','o','markerfacecolor','r','markersize',4','markeredgecolor','r')
hlabels=terlabel('high','medium','low');


%% model comparison

clear all
expnumber = 2;
modVec = [3 1 2 4];
nModels = length(modVec);
modcompidx = 4;
fixedrisk = 0;
MCM = 'BIC';

filename = ['exp' num2str(expnumber) '_cleandata.mat'];
load(filename,'nTrials')
nTrials = sum(nTrials,2);
nSubj = length(nTrials);

if (fixedrisk)
    filepath = ['fits/exp' num2str(expnumber) '_fixedrisk/'];
else
    filepath = ['fits/exp' num2str(expnumber) '/'];
end

nLLMat = nan(nModels,nSubj);
for imodel = 1:nModels
    modidx = modVec(imodel);
    
    load([filepath 'fits_model' num2str(modidx) '.mat'])
    nLLMat(imodel,:) = nLLVec;
    nParamVec(imodel) = size(ML_parameters,2);
    %     AIC.(['model' num2str(imodel)]) = 2*nLLVec + 2*nParamVec(imodel);
end

[AIC, BIC, AICc]= modcomp(nLLMat',nParamVec,nTrials);


% labels and index stuff
modlabels = {'max points','flexible','proportional','min error'};
modcomplabel = modlabels{modcompidx};
modlabels = modlabels(modVec(modVec ~= modcompidx));

nSubj = length(nLLVec);
subtractything = 2-expnumber;
switch MCM
    case 'AIC'
        comparison = bsxfun(@minus,AIC,AIC(:,modVec == modcompidx));
    case 'BIC'
        comparison = bsxfun(@minus,BIC,BIC(:,modVec == modcompidx));
    case 'AICc'
        comparison = bsxfun(@minus,AICc,AICc(:,modVec == modcompidx));
end
comparison(:,modVec == modcompidx) = [];
mediancomp = median(comparison);

% bootstrap the confidence intervals
nBoots = 10000;
for imodel = 1:(nModels-1)
    currVec = comparison(:,imodel);
    sampless = currVec(randi(nSubj,nSubj,nBoots));
    medsamps = sort(median(sampless));
    medCI(imodel,1) = medsamps(.025*nBoots);
    medCI(imodel,2) = medsamps(.975*nBoots);
end

figure;
if nModels == 2
    fill([0 nSubj+1 nSubj+1 0],medCI([1 1 2 2]),...
        0.8*ones(1,3),'EdgeColor','none'); hold on;
    plot([0 nSubj+1],[mediancomp mediancomp],'Color',[0.1 0.1 0.1])
    set(gca,'XTick',[],'XTickLabel',modlabels)
else
    for imodel = 1:(nModels-1)
        fill([imodel-0.475 imodel+0.475 imodel+0.475 imodel-0.475],medCI([imodel imodel imodel+nModels-1 imodel+nModels-1]),...
            0.7*ones(1,3),'EdgeColor','none'); hold on;
        plot([imodel-0.475 imodel+0.475],mediancomp(imodel)*ones(1,2),'Color',[0.1 0.1 0.1])
        set(gca,'XTick',[],'XTick',1:3,'XTickLabel',modlabels)
    end
end

bar(comparison','k')

defaultplot

ylabel(['\Delta ' MCM ' (favoring ' modcomplabel ' model)'])



%% % % % % % % % % % % % % % % % % % % % % %
%       PARAMETER/MODEL RECOVERY
% % % % % % % % % % % % % % % % % % % % % %

% double checking nLLs
% parameter recovery plot


%% simulate data

clear all
expnumber = 2;
imodel = 4;
nSubj = 10;

[simtheta,simdata] = simulate_data(imodel,expnumber,nSubj);


filepath = ['fits/exp' num2str(expnumber) '/'];
save([filepath 'simdata_model' num2str(imodel) '.mat'],'simdata','simtheta')

%% check that none of the subjects have nans

for isubj = 1:nSubj
    cellfun(@(x) sum(isnan(x(:))),simdata{isubj},'UniformOutput',false)
    cellfun(@(x) max(x(:)),simdata{isubj},'UniformOutput',false)
end

%% see what runlist idxs you need to do still
clear all; clc

expnumber = 2;
testmodel = 1;
truemodel = 3;
nSubj = 10;

filepath = ['fits/exp' num2str(expnumber) '/'];

for isubj = 1:nSubj
    load([filepath 'modelrecov_truemodel' num2str(truemodel) '_testmodel' num2str(testmodel) '_subj' num2str(isubj) '.mat'])

    blah = 1:50;
    blah(runlist_completed) = [];
    fprintf([num2str(isubj) ' ' num2str(blah) '\r'])
end


%% look at simulated data

clear all

expnumber = 2;
imodel = 4;

filepath = ['fits/exp' num2str(expnumber) '/'];
load([filepath 'simdata_model' num2str(imodel) '.mat'])

for isubj = 1:10
    xlims = linspace(0,10,11);
    for ipriority = 1:3
        saccerror = hist(simdata{isubj}{ipriority}(:,1),xlims);
        error{ipriority}(isubj,:) = saccerror./sum(saccerror);
        
        if expnumber == 2
            disksize = hist(simdata{isubj}{ipriority}(:,2),xlims);
            discsize{ipriority}(isubj,:) = disksize./sum(disksize);
        end
        
    end
    
    if expnumber == 2
        subplot(1,2,2)
        plot(xlims,disksize,'k-')
        defaultplot;
        xlabel('disc size')
        
        subplot(1,2,1)
    end
    plot(xlims,saccerror,'k-')
    defaultplot
    xlabel('saccade error')
    title(['subj ' num2str(isubj)])
    pause;
end


%% double chek NLL is good for parameter/model recovery

% clear all

truemodel = 4;
expnumber = 2;


testmodelVec = [1 2 3 4];
nModels = length(testmodelVec);

% load true model stuff
filepath = ['fits/exp' num2str(expnumber) '/'];
load([filepath 'simdata_model' num2str(truemodel) '.mat'])

logflag = loadconstraints(truemodel,expnumber);
simtheta(:,logflag) = log(simtheta(:,logflag));

nSubj = 10;
nLLCell = cell(1,nModels+1);
% calculate true nLL
for isubj = 1:nSubj
    nLLCell{1}(isubj) = calc_nLL(truemodel,simtheta(isubj,:),simdata{isubj});
end

for itestmodel = 1:nModels
    itestmodel
    testmodel = testmodelVec(itestmodel);
    
    % load relevant dataset
    load([filepath 'modelrecov_truemodel' num2str(truemodel) '_testmodel' num2str(testmodel) '.mat'])
    
    logflag = loadconstraints(testmodel,expnumber);
    ML_parameters(:,logflag) = log(ML_parameters(:,logflag));
    

    for isubj = 1:nSubj
        % calculate nLL
        nLLCell{itestmodel+1}(isubj) = calc_nLL(testmodel,ML_parameters(isubj,:),simdata{isubj});
    end
    
    
end

% [nLLCell{1}; nLLCell{2}]
[nLLCell{1}; nLLCell{2}; nLLCell{3}; nLLCell{4}; nLLCell{5}]

%% model recovery

clear all

% things to change
expnumber = 1;
modelVec = [3 4];

% things not to change
nSubj = 10;
filepath = ['fits/exp' num2str(expnumber) '/'];
nModels = length(modelVec);
nTrials = sum([250 120 70]);

nLLMat = cell(1,nModels);
nParamMat = cell(1,nModels);
for itruemodel = 1:nModels
    truemodel = modelVec(itruemodel);
    
    nLLMat{itruemodel} = nan(nModels,nSubj);
    nParamMat{itruemodel} = nan(nModels,1);
    for itestmodel = 1:nModels
        testmodel = modelVec(itestmodel);
        
        filename = [filepath 'modelrecov_truemodel' num2str(truemodel) '_testmodel' num2str(testmodel) '.mat'];
        load(filename)
        
        nLLMat{itruemodel}(itestmodel,:) = nLLVec;
        nParamMat{itruemodel}(itestmodel) = size(ML_parameters,2);
        
    end
end

AICMat = cellfun(@(x,y) bsxfun(@plus,2*x,2*y),nLLMat,nParamMat,'UniformOutput',false);
AICcMat = cellfun(@(x,y) bsxfun(@plus,x,(2.*y.*(y+1))./(nTrials-y-1)),AICMat,nParamMat,'UniformOutput',false);
BICMat = cellfun(@(x,y) bsxfun(@plus,2*x,y.*(log(nTrials) - log(2*pi))),nLLMat,nParamMat,'UniformOutput',false);

% confusion matrices
[~,I] = cellfun(@(x) min(x),AICMat,'UniformOutput',false);
AICconfusionMat = nan(nModels);
for itruemodel = 1:nModels
    for itestmodel = 1:nModels
        AICconfusionMat(itruemodel,itestmodel) = sum(I{itruemodel} == itestmodel);
    end
end
AICconfusionMat

[~,I] = cellfun(@(x) min(x),AICcMat,'UniformOutput',false);
AICcconfusionMat = nan(nModels);
for itruemodel = 1:nModels
    for itestmodel = 1:nModels
        AICcconfusionMat(itruemodel,itestmodel) = sum(I{itruemodel} == itestmodel);
    end
end
AICcconfusionMat

[~,I] = cellfun(@(x) min(x),BICMat,'UniformOutput',false);
BICconfusionMat = nan(nModels);
for itruemodel = 1:nModels
    for itestmodel = 1:nModels
        BICconfusionMat(itruemodel,itestmodel) = sum(I{itruemodel} == itestmodel);
    end
end
BICconfusionMat


%% parameter recovery plot

clear all
expnumber = 1;
imodel = 2;
filepath = ['fits/exp' num2str(expnumber) '/'];

load([filepath 'modelrecov_truemodel' num2str(imodel) '_testmodel' num2str(imodel) '.mat'])
bfp = ML_parameters;

load([filepath 'simdata_model' num2str(imodel) '.mat'])
nParams = size(bfp,2);
nSubj = size(bfp,1);

% put things back in log space if they were optimized in log space
[logflag] = loadconstraints(imodel,expnumber);
bfp(:,logflag) = log(bfp(:,logflag));
simtheta(:,logflag) = log(simtheta(:,logflag));

figure;
for iparam = 1:nParams
    subplot(2,3,iparam); hold on;
    for isubj = 1:nSubj;
        plot(bfp(isubj,iparam),simtheta(isubj,iparam),'o'); hold on;
    end
    
    plot([min([bfp(:,iparam);simtheta(:,iparam)]),max([bfp(:,iparam);simtheta(:,iparam)])],...
        [min([bfp(:,iparam);simtheta(:,iparam)]),max([bfp(:,iparam);simtheta(:,iparam)])],'k-')
    xlabel('estimated'); ylabel('actual')
    defaultplot
end

% %% parameter recovery: relationshtip between Jbar and tau
% 
% clear all
% expnumber = 1;
% imodel = 4;
% filepath = ['fits/exp' num2str(expnumber) '/'];
% 
% load([filepath 'modelrecov_truemodel' num2str(imodel) '_testmodel' num2str(imodel) '.mat'])
% bfp = ML_parameters;
% 
% load([filepath 'simdata_model' num2str(imodel) '.mat'])
% nParams = size(bfp,2);
% nSubj = size(bfp,1);
% 
% figure;hold on
% plot([simtheta(:,1) bfp(:,1)]',[simtheta(:,2) bfp(:,2)]','-o');
% 
% figure;
% plot(simtheta(:,1)./simtheta(:,2),bfp(:,1)./bfp(:,2),'o')
% axis([0 50 0 50])


%% make nLL landscape

clear all
imodel = 1;
isubj = 5;

load(['simdata_model' num2str(imodel) '.mat'])
load(['paramrecov_model' num2str(imodel) '.mat'])

thetas = sort([simtheta(isubj,:); ML_parameters(isubj,:)]);

JbartotalVec = linspace(log(thetas(1,1).*0.5),log(thetas(2,1)*2),11);
tauVec = linspace(max(log([ 1e-5 thetas(1,2).*0.5])),max(log([1e-5 thetas(2,2)*2])),11);
% betaVec = linspace(thetas(3).*0.5,thetas(3)*1.5,11);

otherparams = simtheta(isubj,3:end);
for iJ = 1:11
    Jbartotal = JbartotalVec(iJ)
    
    for itau = 1:11
        tau = tauVec(itau);
        
        nLLMat(iJ,itau) = calc_nLL(imodel,[Jbartotal tau otherparams],simdata{isubj});
    end
end

figure;
imagesc(tauVec,JbartotalVec,nLLMat);
hold on;
plot(log(simtheta(isubj,2)),log(simtheta(isubj,1)),'r*')
plot(log(ML_parameters(isubj,2)),log(ML_parameters(isubj,1)),'go')
defaultplot
xlabel('tau'); ylabel('Jbar_{total}')

% logflag = 1:2;
% plb = [0.5 0.01 0.5];
% pub = [20 5 1.5];
% if model == 2
%     plb = [plb 0.3 0];
%     pub = [pub 0.7 0.3];
% end
% lb(logflag) = log(lb(logflag));
% ub(logflag) = log(ub(logflag));
% plb(logflag) = log(plb(logflag));
% pub(logflag) = log(pub(logflag));







%% % % % % % % % % % % % % % % % % % % % % % %
%      MODEL PREDICTION PLOT
% % % % % % % % % % % % % % % % % % % % % % %

% for real data and probability (given ML parameters)
clear all

% ========= simulating a bunch of data per subject =========

expnumber = 1;
imodel = 4;
fixedrisk = [];%'_fixedrisk';
loadpreddata = 0;
indvlplot = 0;

nPriorities = 3;
nTrials = 1e3*ones(1,3); % how many trials to simulate per priority
load(['exp' num2str(expnumber) '_cleandata.mat'],'data')
if (expnumber == 1)
    nSubj = 14;
else
    nSubj = 11;
end
filepath = ['fits/exp' num2str(expnumber) fixedrisk '/'];

% get ML parameter estimate for isubj
load([filepath 'fits_model' num2str(imodel) '.mat'])

if (loadpreddata)
    load([filepath 'modelpred_exp' num2str(expnumber) '_model' num2str(imodel) fixedrisk '.mat'])
else

    preddata = simulate_data(imodel,expnumber,ML_parameters,nTrials);

    
    try
        save([filepath 'modelpred_exp' num2str(expnumber) '_model' num2str(imodel) fixedrisk '.mat'],'preddata','pMat')
    catch
        save([filepath 'modelpred_exp' num2str(expnumber) '_model' num2str(imodel) fixedrisk '.mat'],'preddata')
    end
end

% histograms per subjects
xlims = linspace(0,10,16);

figure;
for isubj = 1:nSubj
    if (indvlplot); figure; end;
    
%     subplot(3,4,isubj);
%     histogram(preddata{isubj}{3}(:,1));
%     title(['subject ' num2str(isubj)])
%     defaultplot
%     xlim([0 15])
    
    for ipriority = 1:nPriorities
        
        % histogram of euclidean error
        datacounts = hist(data{isubj}{ipriority}(:,1),xlims);
        simdatacounts = hist(preddata{isubj}{ipriority}(:,1),xlims);
        error{ipriority}(isubj,:) = datacounts./sum(datacounts);
        simerror{ipriority}(isubj,:) = simdatacounts./sum(simdatacounts);
        
        if (indvlplot)
            subplot(3,2,2*ipriority-1)
            plot(xlims,error{ipriority}(isubj,:),'k')
            hold on;
            plot(xlims,simerror{ipriority}(isubj,:),'Color',aspencolors('booger'));
            defaultplot
            if ipriority == 1; title('euclidean error'); end
        end
        
        if (expnumber == 2)
            % histogram of disc size
            datacounts = hist(data{isubj}{ipriority}(:,2),xlims);
            simdatacounts = hist(preddata{isubj}{ipriority}(:,2),xlims);
            discsize{ipriority}(isubj,:) = datacounts./sum(datacounts);
            simdiscsize{ipriority}(isubj,:)  = simdatacounts./sum(simdatacounts);
            
            if (indvlplot)
                subplot(3,2,2*ipriority)
                plot(xlims,discsize{ipriority}(isubj,:),'k')
                hold on;
                plot(xlims,simdiscsize{ipriority}(isubj,:),'Color',aspencolors('booger'));
                defaultplot
                if ipriority == 1; title('disc size'); end
            end
        end
    end
    %     pause;
end

% =========== group plot =====================

figure;
colorMat = {'r','b','k'};
if (expnumber == 2)
    ha = tight_subplot(3,3,{[.03 .03],[.03 .07]},[.1 .01],[.1 .01]);
else
    ha = tight_subplot(1,3,.03,[.26 .05],[.11 .05]);
end

meanerror = cellfun(@mean,error,'UniformOutput',false);
semerror = cellfun(@(x) std(x)./sqrt(size(x,1)),error,'UniformOutput',false);
meansimerror = cellfun(@mean,simerror,'UniformOutput',false);
semsimerror = cellfun(@(x) std(x)./sqrt(size(x,1)),simerror,'UniformOutput',false);

if (expnumber == 2)
    meandiscsize = cellfun(@mean,discsize,'UniformOutput',false);
    semdiscsize = cellfun(@(x) std(x)./sqrt(size(x,1)),discsize,'UniformOutput',false);
    meansimdiscsize = cellfun(@mean,simdiscsize,'UniformOutput',false);
    semsimdiscsize = cellfun(@(x) std(x)./sqrt(size(x,1)),simdiscsize,'UniformOutput',false);
end

for ipriority = 1:nPriorities
    
    if(expnumber == 2)
        axes(ha(3*ipriority-2))
    else
        axes(ha(ipriority))
    end
    fill([xlims fliplr(xlims)],[meansimerror{ipriority}-semsimerror{ipriority}...
        fliplr(meansimerror{ipriority}+semsimerror{ipriority})],colorMat{ipriority},'EdgeColor','none','FaceAlpha',0.4);
    hold on;
    errorbar(xlims,meanerror{ipriority},semerror{ipriority},'Color','k','LineStyle','none','LineWidth',1);
    defaultplot
    axis([0 10 0 0.4])
    
    if expnumber == 2
        %          axis([0 10 0 0.6])
        if ipriority == 3
            xlabel('error');
        else
            set(ha(3*ipriority-2),'XTickLabel','');
        end
        ylabel('proportion','FontSize',14);
        
    else
        
        xlabel('error','FontSize',16)
        set(ha(ipriority),'YTick',[0 0.2 0.4],'FontSize',12);
        if ipriority ~= 1
            set(ha(ipriority),'YTickLabel','');
        else
            ylabel('proportion','FontSize',16)
        end
        
    end
    
    
    if (expnumber == 2)
        % discsize
        axes(ha(3*ipriority-1))
        fill([xlims fliplr(xlims)],[meansimdiscsize{ipriority}-semsimdiscsize{ipriority}...
            fliplr(meansimdiscsize{ipriority}+semsimdiscsize{ipriority})],colorMat{ipriority},'EdgeColor','none','FaceAlpha',0.4);
        hold on;
        errorbar(xlims,meandiscsize{ipriority},semdiscsize{ipriority},'Color','k','LineStyle','none','LineWidth',1);
        defaultplot
        axis([0 10 0 0.6])
        if ipriority == 3
            xlabel('disc size','FontSize',14);
        else
            set(ha(3*ipriority-1),'XTickLabel','');
        end
        set(ha(3*ipriority-1),'YTickLabel','');
        
    end
    
end

% ========== quantile correlation plot per subject ===========

nQuants = 6;
for isubj = 1:11
    if (indvlplot); figure; end
    
    for ipriority = 1:nPriorities
        currdata = data{isubj}{ipriority}(:,1);
        [currdata,idx] = sort(currdata);
        quantVec = round(linspace(0,length(currdata),nQuants+1));
        
        currsimdata = preddata{isubj}{ipriority}(:,1);
        [currsimdata,simidx] = sort(currsimdata);
        simquantVec = round(linspace(0,length(currsimdata),nQuants+1));
        for iquant = 1:nQuants
            meanquanterror{ipriority}(isubj,iquant) = mean(currdata(quantVec(iquant)+1:quantVec(iquant+1)));
            meanquantdiscsize{ipriority}(isubj,iquant) = mean(data{isubj}{ipriority}(idx(quantVec(iquant)+1:quantVec(iquant+1)),2));
            
            meanquantsimerror{ipriority}(isubj,iquant) = mean(currsimdata(simquantVec(iquant)+1:simquantVec(iquant+1)));
            meanquantsimdiscsize{ipriority}(isubj,iquant) = mean(preddata{isubj}{ipriority}(simidx(simquantVec(iquant)+1:simquantVec(iquant+1)),2));
        end
        
        if (indvlplot)
            subplot(3,1,ipriority)
            plot(meanquanterror{ipriority}(isubj,:),meanquantdiscsize{ipriority}(isubj,:),'k');
            hold on;
            plot(meanquantsimerror{ipriority}(isubj,:),meanquantsimdiscsize{ipriority}(isubj,:),'Color',aspencolors('booger'));
        end
    end
    if (indvlplot); pause; end
    
end


% ================ group plot ====================

meanmeanquanterror = cellfun(@mean,meanquanterror,'UniformOutput',false);
meanmeanquantdiscsize = cellfun(@mean,meanquantdiscsize,'UniformOutput',false);
semmeanquantdiscsize= cellfun(@(x) std(x)./sqrt(size(x,1)),meanquantdiscsize,'UniformOutput',false);
meanmeanquantsimerror = cellfun(@nanmean,meanquantsimerror,'UniformOutput',false);
meanmeanquantsimdiscsize = cellfun(@mean,meanquantsimdiscsize,'UniformOutput',false);
semmeanquantsimdiscsize= cellfun(@(x) std(x)./sqrt(size(x,1)),meanquantsimdiscsize,'UniformOutput',false);

% figure;
colorMat = {'r','b','k'};
for ipriority = 1:nPriorities
    axes(ha(3*ipriority))
    %         subplot(3,3,6+ipriority)
    hold on
    plot_summaryfit(meanmeanquantsimerror{ipriority},[],[],meanmeanquantsimdiscsize{ipriority},...
        semmeanquantsimdiscsize{ipriority},[],colorMat{ipriority})
    plot_summaryfit(meanmeanquanterror{ipriority},meanmeanquantdiscsize{ipriority},semmeanquantdiscsize{ipriority},...
        [],[],'k')
    
    ylabel('disc size');
    axis([0 6 2 4])
    set(ha(3*ipriority),'XTick',[0 3 6],'YTick',[2 3 4]);
    set(ha(3*ipriority),'YTickLabel',[2 3 4])
    if (expnumber == 1)
        set(ha(ipriority),'XTickLabel',[0 3 6])
        xlabel('error');
    end
end
if (expnumber == 2);set(ha(3*ipriority),'XTickLabel',[0 3 6]); end
xlabel('error');







%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%           RANDOM THINGS
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear all
expnumber = 2;
model = 2;

load(['exp' num2str(expnumber) '_cleandata.mat'])
filepath = ['fits/exp' num2str(expnumber) '/'];
load([filepath 'fits_model' num2str(model) '.mat'])


%% look at LL as a function of one parameter

isubj = 1;


bfp = ML_parameters(isubj,:);
propVec = [.5 .7 .9 1 1.1 1.3 1.5];
nProps = length(propVec);

figure;
for iparam = 1:6
    iparam
    
    nLL = nan(1,nProps);
    for iprop = 1:nProps
        
        testtheta = bfp;
        testtheta(iparam) = bfp(iparam)*propVec(iprop);
        testtheta(1:2) = log(testtheta(1:2));
        
        nLL(iprop) = calc_nLL(model,testtheta,data{isubj});
    end
    
    subplot(2,3,iparam)
    plot(propVec,nLL,'k-');
    defaultplot;
    title(num2str(iparam))
end


%% looking at LL calculation as a function of number of grids for gamma distribution

iparam = 1;
isubj = 1;

gridVec = linspace(100,1000,100);
nGrids = length(gridVec);

MLtheta = 0;
if (MLtheta)
    theta = ML_parameters(isubj,:);
    propVec = [1];
    nProps = length(propVec);
else
    MU = mean(ML_parameters);
    SIGMA = cov(ML_parameters);
    theta = abs(mvnrnd(MU,SIGMA));
end

tic;
nLLMat = nan(nProps,nGrids);
for iprop = 1:nProps
    
    testtheta = theta;
    testtheta(iparam) = theta(iparam)*propVec(iprop);
    testtheta(1:2) = log(testtheta(1:2));
    
    for igrid = 1:nGrids
        if ~mod(igrid,10)
            disp(nLLMat(iprop,igrid-1))
        end
        
        nLLMat(iprop,igrid) = calc_nLL(model,testtheta,data{isubj},gridVec(igrid));
    end
end
toc;
figure;
plot(gridVec,nLLMat);
defaultplot
xlabel('number of grids'); ylabel('negative log-likelihood')
if (MLtheta)
    title(['ML \theta for subject ' num2str(isubj)])
else
    title(['\theta = ' num2str(theta)])
end
% imagesc(nLLMat)

% dim = [0.2 0.5 0.3 0.3];
% str = {'Straight Line Plot','from 1 to 10'};
% annotation('textbox',dim,'String',str,'FitBoxToText','on');

%% looking at LL calculation as a function of number of grids for rVec


isubj = 1;

gridVec = linspace(100,1000,100);
nGrids = length(gridVec);

MLtheta = 1;
if (MLtheta)
    theta = ML_parameters(isubj,:);
else
    MU = mean(ML_parameters);
    SIGMA = cov(ML_parameters);
    theta = abs(mvnrnd(MU,SIGMA));
end

nLLMat = nan(1,nGrids);
testtheta = theta;
testtheta(1:2) = log(testtheta(1:2));
for igrid = 1:nGrids
    if ~mod(igrid,10)
        disp(nLLMat(igrid-1))
    end
    
    nLLMat(igrid) = calc_nLL(model,testtheta,data{isubj},gridVec(igrid));
end

figure;
plot(gridVec,nLLMat);
defaultplot
xlabel('number of grids'); ylabel('negative log-likelihood')
if (MLtheta)
    title(['ML \theta for subject ' num2str(isubj)])
else
    title(['\theta = ' num2str(theta)])
end



%% nLL as a function of numerical integration samples
sampVec = [100 200 500 1000 10000 15000];
isubj = 4;
for isamp = 1:(length(sampVec)-1)
    isamp
    nLLVec(isamp) = calc_nLL(testmodel,[log(ML_parameters(isubj,1:2)) ML_parameters(isubj,3:end)],simdata{isubj},sampVec(isamp));
end
figure;
plot(log(sampVec),nLLVec)

%% plot bias in polar angle as a function of location
% exp 1

% load group_data.mat from other SWM folder

clear all
load('/Users/blobface/Research/SWM/VWM-Eye-Movements-Uncertainty/group_data.mat')

% response
r_x = group_data(:,6);
r_y = group_data(:,7);
r_angle = atand(r_x./r_y); % response in polar angles

idx = (r_x < 0) & (r_y < 0); % third quadrant, both negative
r_angle(idx) = r_angle(idx) + 180;
idx = (r_x < 0) & (r_y > 0); % second quadrant, negative x
r_angle(idx) = r_angle(idx) + 180;
idx = (r_x > 0) & (r_y < 0); % fourth quadrant, negative y
r_angle(idx) = r_angle(idx) + 360;

% target
t_x = group_data(:,8);
t_y = group_data(:,9);
t_angle = atand(t_x./t_y);

idx = (t_x < 0) & (t_y < 0); % third quadrant, both negative
t_angle(idx) = t_angle(idx) + 180;
idx = (t_x < 0) & (t_y > 0); % second quadrant, negative x
t_angle(idx) = t_angle(idx) + 180;
idx = (t_x > 0) & (t_y < 0); % fourth quadrant, negative y
t_angle(idx) = t_angle(idx) + 360;

idx = isnan(t_angle);
idx = idx + isnan(r_angle);
idx = logical(idx);

t_angle(idx) = [];
r_angle(idx) = [];

t_angle_rad = deg2rad(t_angle);
r_angle_rad = deg2rad(r_angle);

biass = circ_dist(r_angle_rad, t_angle_rad);

figure; plot(t_angle,'o')

figure; plot(t_angle,biass,'o')

%% average and plot

angless = unique(round(t_angle));
nAngles = length(angless);

for iangle = 1:nAngles
    anglee = angless(iangle);
    
    idx = round(t_angle) == anglee;
    
    biasVec(iangle) = mean(biass(idx));
    
end

figure
plot(angless,biasVec,'o')

%% pick every 10 angles

angless = [10:10:80 100:10:170 180:10:260 280:10:350];
nAngles = length(angless);

biasVec = [];
for iangle = 1:nAngles
    anglee = angless(iangle);
    
    idx = round(t_angle) == anglee;
    
    biasVec(iangle) = mean(biass(idx));
    
end

figure
plot(angless,biasVec,'o')

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%       CLUSTER RELATED
% % % % % % % % % % % % % % % % % % % % % % % % % % % 

clear all
expnumber = 2;
modelnum = 1;
if expnumber == 1; nSubj = 14; else; nSubj = 11;end

filepath = ['fits/exp' num2str(expnumber) '/'];

for isubj = 1:nSubj
    load([filepath 'fits_model' num2str(modelnum) '_subj' num2str(isubj) '.mat'],'runlist_completed')
    runlist.(['subj' num2str(isubj)]) = sort(runlist_completed);
end