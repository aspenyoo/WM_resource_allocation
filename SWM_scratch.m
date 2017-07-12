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

%% % % % % % % % % % % % % % % % % % % %
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



%% % % % % % % % % % % % % % % % % % % % % % % % % %
%       MODEL RELATED
% % % % % % % % % % % % % % % % % % % % % % % % % %

%% understanding Jbar, tau, and sigma

Jbar = 1;
tau = .1;
beta = 0.01;

[JVec] = loadvar({'JVec',Jbar,tau});
% JVec = linspace(1,2,100);
Jpdf = gampdf(JVec,Jbar/tau,tau);
Jpdf = Jpdf./qtrapz(Jpdf); % normalize


% figure;
subplot(3,1,1)
plot(JVec,Jpdf,'k-'); defaultplot;hold on
subplot(3,1,2)
plot(JVec./Jbar,Jpdf,'k-'); defaultplot;hold on
subplot(3,1,3)
plot(1./sqrt(JVec),Jpdf,'k-'); defaultplot; hold on

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
beta = 1;
N = 40;
JbartotalVec = linspace(1e-3,20,N);

pVec = nan(N,3);
for iJbartotal = 1:N
    iJbartotal
    
    Jbartotal = JbartotalVec(iJbartotal);
    pVec(iJbartotal,:) = calc_optimal_pVec([Jbartotal tau beta]);
end

plot(bsxfun(@times,pVec,JbartotalVec'))

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




%% % % % % % % % % % % % % % % % % % % % % % % %
%     GET ML PARAMETER ESTIMATES
% % % % % % % % % % % % % % % % % % % % % % % %
clear all

imodel = 3;
testmodel = 3;
fakedata = 0;
expnumber = 2;
isriskfixed = 0;

subjVec = 1:11;
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

% for isubj = 1:nSubj;
    isubj = 14;
    
%     load(['fits/exp' num2str(expnumber) '/fits_model' num2str(testmodel) '_subj' num2str(isubj) '.mat'])
    runlist = 1:runmax;
%     runlist(unique(runlist_completed)) = [];
    
    fit_parameters(testmodel,isubj,runlist,runmax,truemodel,expnumber)
% end

%% % % % % % % % % % % % % % % % % % % % % % %
%       REAL DATA STUFF
% % % % % % % % % % % % % % % % % % % % % % %



%% check nLL

clear all
expnumber = 2;
imodel = 1;
subjVec = [1 3:11];

nSubj = 10;
load('cleandata.mat')

filepath = ['fits/exp' num2str(expnumber) '/'];
load([filepath 'fits_model' num2str(imodel) '.mat'])
ML_parameters(:,1:2) = log(ML_parameters(:,1:2));

nLL2 = nan(1,nSubj);
for isubj = 1:nSubj
    isubj
    
    nLL2(isubj) = calc_nLL(imodel,ML_parameters(isubj,:),data{subjVec(isubj)});
end

[nLLVec; nLL; nLL2]

%% optimal pVec for model 1

clear all
imodel = 1;

load(['fits_model' num2str(imodel) '.mat'])
nSubj = 11;

pMat = nan(nSubj,3);
for isubj = 1:nSubj
    isubj
    pMat(isubj,:) = calc_optimal_pVec(ML_parameters(isubj,:));
end

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
fixedrisk = '_fixedrisk';
filepath = ['fits/exp' num2str(expnumber)  fixedrisk '/'];
load([filepath 'fits_model' num2str(imodel) '.mat'])
nSubj = size(ML_parameters,1);

pMat = ML_parameters(:,end-1:end);
pMat(:,3) = 1-sum(pMat,2);

% plot(pMat)
figure
plot(bsxfun(@times,pMat,ML_parameters(:,1)))
hold on;

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

%% nTrials for each subject
clear all
expnumber = 1;

switch expnumber
    case 1
        filename = 'cleandata_nodisc.mat';
    case 2
        filename = 'cleandata.mat';
end
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

%% model comparison

clear all
expnumber = 2;
nModels = expnumber+1;
modcompidx = 2;
fixedrisk = 0;
MCM = 'BIC';

switch expnumber
    case 1
        filename = 'cleandata_nodisc.mat';
    case 2
        filename = 'cleandata.mat';
end
load(filename,'nTrials')
nTrials = sum(nTrials,2);

if (fixedrisk)
    filepath = ['fits/exp' num2str(expnumber) '_fixedrisk/'];
else
    filepath = ['fits/exp' num2str(expnumber) '/'];
end

for imodel = (4-nModels):3
    load([filepath 'fits_model' num2str(imodel) '.mat'])
    nLL.(['model' num2str(imodel)]) = nLLVec;
    nParamVec(imodel) = size(ML_parameters,2);
    %     AIC.(['model' num2str(imodel)]) = 2*nLLVec + 2*nParamVec(imodel);
end
if expnumber == 1; nParamVec(1) = []; end

nLLMat = reshape(struct2mat(nLL),[length(nTrials),nModels])';
[AIC  BIC  AICc ]= modcomp(nLLMat',nParamVec,nTrials);


% labels and index stuff
modlabels = {'optimal','free','fixed'};
modcomplabel = modlabels{modcompidx};
modidx = (4-nModels):3;
modidx(modidx == modcompidx) = [];
modlabels = modlabels(modidx);

nSubj = length(nLLVec);
switch MCM
    case 'AIC'
        comparison = bsxfun(@minus,AIC,AIC(:,modcompidx));
    case 'BIC'
        comparison = bsxfun(@minus,BIC,BIC(:,modcompidx));
    case 'AICc'
        comparison = bsxfun(@minus,AICc,AICc(:,modcompidx));
end
comparison(:,modcompidx) = [];
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
if nModels == 3
    for imodel = 1:(nModels-1)
        %         imodel = modidx(imodel);
        fill([imodel-0.475 imodel+0.475 imodel+0.475 imodel-0.475],medCI([1 1 2 2]),...
            0.7*ones(1,3),'EdgeColor','none'); hold on;
        plot([imodel-0.475 imodel+0.475],mediancomp(imodel)*ones(1,2),'Color',[0.1 0.1 0.1])
        set(gca,'XTick',[],'XTick',[1 2],'XTickLabel',modlabels)
    end
else
    fill([0 nSubj+1 nSubj+1 0],medCI([1 1 2 2]),...
        0.8*ones(1,3),'EdgeColor','none'); hold on;
    plot([0 nSubj+1],[mediancomp mediancomp],'Color',[0.1 0.1 0.1])
    set(gca,'XTick',[],'XTickLabel',modlabels)
end

bar(comparison','k')

defaultplot

ylabel(['\Delta AIC (favoring ' modcomplabel ' model)'])



%% % % % % % % % % % % % % % % % % % % % % %
%       PARAMETER/MODEL RECOVERY
% % % % % % % % % % % % % % % % % % % % % %

% double checking nLLs
% parameter recovery plot


%% simulate data

clear all
expnumber = 2;
imodel = 3;
nSubj = 10;

% switch model
%     case 1
%         logflag = logical([1 1 0]);
%     case 2
%         logflag = logical([1 1 0 0 0]);
% end

filepath = ['fits/exp' num2str(expnumber) '/'];
load([filepath 'fits_model' num2str(imodel) '.mat'])
% ML_parameters(logflag) = log(ML_parameters(logflag));
MU = mean(ML_parameters);
SIGMA = cov(ML_parameters);

simtheta = mvnrnd(MU,SIGMA,nSubj);
simtheta = abs(simtheta); % hacky way to enforce positive parameter values

% simtheta(:,logflag) = exp(simtheta(:,logflag));

nTrials = [250 120 70]; % mean number of trials across actual participants
for isubj = 1:nSubj
    isubj
    simdata{isubj} = simulate_data(imodel,expnumber, simtheta(isubj,:),nTrials);
end
save([filepath 'simdata_model' num2str(imodel) '.mat'],'simdata','simtheta')


%% look at simulated data

clear all

expnumber = 2;
imodel = 2;

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

clear all

truemodel = 3;
expnumber = 2;

testmodelVec = [3];
nModels = length(testmodelVec);

% load true model stuff
filepath = ['fits/exp' num2str(expnumber) '/'];
load([filepath 'simdata_model' num2str(truemodel) '.mat'])


nSubj = 10;
nLLCell = cell(1,nModels);
% calculate true nLL
for isubj = 1:nSubj
    nLLCell{1}(isubj) = calc_nLL(truemodel,[log(simtheta(isubj,1:2)) simtheta(isubj,3:end)],simdata{isubj});
end

for itestmodel = 1:nModels
    itestmodel
    testmodel = testmodelVec(itestmodel);
    
    % load relevant dataset
    load([filepath 'modelrecov_truemodel' num2str(truemodel) '_testmodel' num2str(testmodel) '.mat'])
    for isubj = 1:nSubj
        % calculate nLL
        nLLCell{itestmodel+1}(isubj) = calc_nLL(testmodel,[log(ML_parameters(isubj,1:2)) ML_parameters(isubj,3:end)],simdata{isubj});
    end
    
    
end

[nLLCell{1}; nLLCell{2}]

%% model recovery

clear all

% things to change
expnumber = 2;
modelVec = [2 3];

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

% which one wins
[M,I] = cellfun(@(x) min(x),AICMat,'UniformOutput',false);
sum(I{1} == 1)
sum(I{2} == 2)

[M,I] = cellfun(@(x) min(x),AICcMat,'UniformOutput',false);
sum(I{1} == 1)
sum(I{2} == 2)

[M,I] = cellfun(@(x) min(x),BICMat,'UniformOutput',false);
sum(I{1} == 1)
sum(I{2} == 2)



%% parameter recovery plot

clear all
expnumber = 2;
imodel = 2;
filepath = ['fits/exp' num2str(expnumber) '/'];

load([filepath 'modelrecov_truemodel' num2str(imodel) '_testmodel' num2str(imodel) '.mat'])
bfp = ML_parameters;

load([filepath 'simdata_model' num2str(imodel) '.mat'])
nParams = size(bfp,2);

figure;
for iparam = 1:nParams
    subplot(2,3,iparam);
    plot(bfp(:,iparam),simtheta(:,iparam),'ko'); hold on;
    plot([min([bfp(:,iparam);simtheta(:,iparam)]),max([bfp(:,iparam);simtheta(:,iparam)])],...
        [min([bfp(:,iparam);simtheta(:,iparam)]),max([bfp(:,iparam);simtheta(:,iparam)])],'k-')
    xlabel('estimated'); ylabel('actual')
    defaultplot
end

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
imodel = 2;
fixedrisk = [];%'_fixedrisk';
loadpreddata = 0;
indvlplot = 0;

nPriorities = 3;
nTrials = 1e3*ones(1,3); % how many trials to simulate per priority
if (expnumber == 1)
    load('cleandata_nodisc.mat','data')
    nSubj = 14;
else
    load('cleandata.mat','data')
    nSubj = 11;
end
filename = ['fits/exp' num2str(expnumber) fixedrisk '/'];

% get ML parameter estimate for isubj
load([filename 'fits_model' num2str(imodel) '.mat'])

if (loadpreddata)
    load([filename 'modelpred_exp' num2str(expnumber) '_model' num2str(imodel) fixedrisk '.mat'],'preddata')
else
    preddata = cell(1,nSubj);
    for isubj = 1:nSubj
        isubj
        Theta = ML_parameters(isubj,:);
        preddata{isubj} = simulate_data(imodel,expnumber,Theta,nTrials);
    end
    
    save([filename 'modelpred_exp' num2str(expnumber) '_model' num2str(imodel) fixedrisk '.mat'],'preddata')
end

% histograms per subjects
xlims = linspace(0,10,16);

for isubj = 1:nSubj
    if (indvlplot); figure; end;
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

if (expnumber == 2)
    load('cleandata.mat')
else
    load('cleandata_nodisc.mat')
end
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


plot(t_angle,biass,'o')

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