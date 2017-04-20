%% understanding Jbar, tau, and sigma

Jbar = 10;
tau = .01;
beta = 0.01;

[JVec] = loadvar({'JVec',Jbar,tau});
% JVec = linspace(1,2,100);
Jpdf = gampdf(JVec,Jbar/tau,tau);
% Jpdf = Jpdf./qtrapz(Jpdf); % normalize

p_r = calc_pdf_r(beta,Jbar);

% figure; 
subplot(2,1,1)
plot(JVec,Jpdf,'k-'); defaultplot;hold on
subplot(2,1,2)
plot(1./sqrt(JVec),Jpdf,'k-'); defaultplot; hold on



%% ==========================================================
%             STUFF WITH ACTUAL DATA!
% ===========================================================
% 02.5.2017

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

%% % % % % % % % % % % % % % % % % % % % % % % %
%     GET ML PARAMETER ESTIMATES 
% % % % % % % % % % % % % % % % % % % % % % % % 

imodel = 2;
nSubj = 10;
fakedata = 1;
expnumber = 1;

filepath = ['fits/exp' num2str(expnumber) '/'];
if (fakedata)
    pretxt = 'paramrecov';
else
    pretxt = 'fits';
end

switch imodel
    case {1,3}
        nParams = 3;
    case 2 
        nParams = 5;
end
if (expnumber == 1); nParams = nParams - 1; end

bfp = nan(nSubj,nParams);
nLL = nan(1,nSubj);
for isubj = 1:nSubj
    isubj;
    load([filepath pretxt '_model' num2str(imodel) '_subj' num2str(isubj) '.mat'])
    blah = ML_parameters(nLLVec == min(nLLVec),:);
    bfp(isubj,:) = blah(1,:);
    nLL(isubj) = min(nLLVec);
end
ML_parameters = bfp;
nLLVec = nLL;
save([filepath pretxt '_model' num2str(imodel) '.mat'],'ML_parameters','nLLVec')

%% parameter recovery plot

clear all
expnumber = 1;
imodel = 3;
filepath = ['fits/exp' num2str(expnumber) '/'];

load([filepath 'paramrecov_model' num2str(imodel) '.mat'])
bfp = ML_parameters;

load([filepath 'simdata_model' num2str(imodel) '.mat'])
nParams = size(bfp,2);

for iparam = 1:nParams
    subplot(2,3,iparam);
    plot(bfp(:,iparam),simtheta(:,iparam),'ko'); hold on; 
    plot([min([bfp(:,iparam);simtheta(:,iparam)]),max([bfp(:,iparam);simtheta(:,iparam)])],...
        [min([bfp(:,iparam);simtheta(:,iparam)]),max([bfp(:,iparam);simtheta(:,iparam)])],'k-')
    defaultplot
end

%% double chek NLL is good for parameter recovery

clear all

imodel = 1;

load(['simdata_model' num2str(imodel) '.mat'])
load(['paramrecov_model' num2str(imodel) '.mat'])
bfp = ML_parameters;
bfp(:,1:2) = log(bfp(:,1:2));

nSubj = 10;
% [nLLVec2, nLLVec3] = deal(nan(1,nSubj));
for isubj = 1:nSubj
    isubj
    
    nLLVec4(isubj) = calc_nLL(imodel,bfp(isubj,:),simdata{isubj});
    nLLVec3(isubj) = calc_nLL(imodel,[log(simtheta(isubj,1:2)) simtheta(isubj,3:end)],simdata{isubj});
end

[nLLVec; nLLVec3; nLLVec4]

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

%% optimal pVec for model 1

clear all
imodel = 1;

load(['fits_model' num2str(imodel) '.mat'])
nSubj = 11;

pMat = nan(nSubj,3);
for isubj = 1:nSubj
    
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
filepath = ['fits/exp' num2str(expnumber) '/'];
load([filepath 'fits_model' num2str(imodel) '.mat'])
nSubj = size(ML_parameters,1);

pMat = ML_parameters(:,end-1:end);
pMat(:,3) = 1-sum(pMat,2);

% plot(pMat)
figure
plot(bsxfun(@times,pMat,ML_parameters(:,1)))
hold on;


%% plot triangle (ternary) plot

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

%% double chek NLL is good

imodel = 1;

load(['fits_model' num2str(imodel) '.mat'])
load('cleandata.mat','data')
bfp = ML_parameters;
bfp(:,1:2) = log(bfp(:,1:2));

nSubj = 11;
nLLVec2 = nan(1,nSubj);
for isubj = 1:nSubj
    isubj
    
    nLLVec2(isubj) = calc_nLL(imodel,bfp(isubj,:),data{isubj});
end

[nLLVec; nLLVec2]

%% model comparison

clear all
expnumber = 1;
filepath = ['fits/exp' num2str(expnumber) '/'];

nParamVec = nan(1,2);
for imodel = 2:3
    load([filepath 'fits_model' num2str(imodel) '.mat'])
    nLL.(['model' num2str(imodel)]) = nLLVec;
    nParamVec(imodel) = size(ML_parameters,2);
    AIC.(['model' num2str(imodel)]) = 2*nLLVec + 2*nParamVec(imodel);
end

modcompidx = 3;
nSubj = length(nLLVec);
comparison = structfun(@(x) x-AIC.(['model' num2str(modcompidx)]),AIC,'UniformOutput',false);
comparison = comparison.model2;
meancomp = mean(comparison);
semcomp = std(comparison)/sqrt(length(comparison));

figure;
fill([0 nSubj+1 nSubj+1 0],[meancomp-semcomp meancomp-semcomp meancomp+semcomp meancomp+semcomp],...
    0.7*ones(1,3),'EdgeColor','none'); 
hold on;
bar(comparison,'k')

defaultplot
set(gca,'XTick',[],'XTickLabel',[])
ylabel(['\Delta AIC (favoring model ' num2str(modcompidx) ')'])


%% % % % % % % % % % % % % % % % % % % % % % %
%       PLOTS 
% % % % % % % % % % % % % % % % % % % % % % % 

% for real data and one simulated dataset (using ML parameters)
% 1. disc size as a function of euclidean error
% 2. histogram of error and disc sizes for each priority
% 3. median and IQR for error and disc size for each priority
 
% for real data and probability (given ML parameters)
% 4. data overlayed onto heatmap of probabilities

%% simulate one dataset

close all

expnumber = 1;
nPriorities = 3;
imodel = 3;
isubj = 1; % fitted data isubj = actual data isubj - 3
if (expnumber == 1)
    load('cleandata_nodisc.mat','data')
else
    load('cleandata.mat','data')
end

% get ML parameter estimate for isubj
filename = ['fits/exp' num2str(expnumber) '/'];
load([filename 'fits_model' num2str(imodel) '.mat'])
Theta = ML_parameters(isubj,:); 

nTrials = nan(1,3);
for ipriority = 1:nPriorities
    nTrials(ipriority) = size(data{isubj}{ipriority},1);
end

simdata = simulate_data(imodel,expnumber,Theta,nTrials);
%% 1. plot disc size as a function of euclidean error

% plot simulated data on top of real data
for ipriority = 1:nPriorities
    subplot(2,2,ipriority)
    plot(data{isubj}{ipriority}(:,1),data{isubj}{ipriority}(:,2),'ko');
    hold on;
    plot(simdata{ipriority}(:,1),simdata{ipriority}(:,2),'r*');
    defaultplot
end

%% 2. plot histogram of error and disc size for each priority

figure;
for ipriority = 1:nPriorities
    
    % plot error
    subplot(3,2,2*ipriority-1)
    hist([data{isubj}{ipriority}(:,1) simdata{ipriority}(:,1)])
    
    if (expnumber == 2)
    % plot discsize
    subplot(3,expnumber,expnumber*ipriority)
    hist([data{isubj}{ipriority}(:,2) simdata{ipriority}(:,2)])
    end
end

%% 3. plot median and IQR for error and disc size for each priority

figure;
for ipriority = 1:nPriorities
    
    % DISTANCE ERROR
    % error from real data
    subplot(1,2,1)
    sorteddata = sort(data{isubj}{ipriority}(:,1));
    quart25 = sorteddata(round(length(sorteddata)*0.25));
    quart75 = sorteddata(round(length(sorteddata)*0.75));   
    med = median(sorteddata);
    hold on; errorbar(ipriority,med,med-quart25,quart75-med,'Marker','.','Color','k');
    
    % simualted data error
    sorteddata = sort(simdata{ipriority}(:,1));
    quart25 = sorteddata(round(length(sorteddata)*0.25));
    quart75 = sorteddata(round(length(sorteddata)*0.75));   
    med = median(sorteddata);
    hold on; errorbar(ipriority+0.2,med,med-quart25,quart75-med,'Marker','.','Color','r');
    
    % DISC SIZE
    % real data
    subplot(1,2,2)
    sorteddata = sort(data{isubj}{ipriority}(:,2));
    quart25 = sorteddata(round(length(sorteddata)*0.25));
    quart75 = sorteddata(round(length(sorteddata)*0.75));   
    med = median(sorteddata);
    hold on; errorbar(ipriority,med,med-quart25,quart75-med,'Marker','.','Color','k');
    
    % simulated disc size
    sorteddata = sort(simdata{ipriority}(:,2));
    quart25 = sorteddata(round(length(sorteddata)*0.25));
    quart75 = sorteddata(round(length(sorteddata)*0.75));   
    med = median(sorteddata);
    hold on; errorbar(ipriority+0.2,med,med-quart25,quart75-med,'Marker','.','Color','r');
end

subplot(1,2,1);
defaultplot
xlim([0.5 3.5])
title('saccadic error (dva)')

subplot(1,2,2)
defaultplot
xlim([0.5 3.5])
title('disc size (dva)')

%% simulating a bunch of data per subject 
% 4. plot disc size as a function of euclidean error
clear all

expnumber = 1;
nPriorities = 3;
imodel = 3;
nTrials = 1e3*ones(1,3); % how many trials to simulate per priority
if (expnumber == 1)
    load('cleandata_nodisc.mat','data')
    nSubj = 14;
else
    load('cleandata.mat','data')
    nSubj = 11;
end
filename = ['fits/exp' num2str(expnumber) '/'];

% get ML parameter estimate for isubj
load([filename 'fits_model' num2str(imodel) '.mat'])

loadpreddata = 1;
if (loadpreddata)
    load([filename 'modelpred_exp' num2str(expnumber) '_model' num2str(imodel) '.mat'],'preddata')
else
    preddata = cell(1,nSubj);
    for isubj = 1:nSubj
        isubj
        Theta = ML_parameters(isubj,:);
        preddata{isubj} = simulate_data(imodel,expnumber,Theta,nTrials);
    end
    
    save([filename 'modelpred_exp' num2str(expnumber) '_model' num2str(imodel) '.mat'],'preddata')
end
%% histograms per subjects
xlims = linspace(0,10,11);

for isubj = 1:nSubj
    figure;
    for ipriority = 1:nPriorities
        
        % histogram of euclidean error
        datacounts = hist(data{isubj}{ipriority}(:,1),xlims);
        simdatacounts = hist(preddata{isubj}{ipriority}(:,1),xlims);
        error{ipriority}(isubj,:) = datacounts./sum(datacounts);
        simerror{ipriority}(isubj,:) = simdatacounts./sum(simdatacounts);
        
        subplot(3,2,2*ipriority-1)
        plot(xlims,error{ipriority}(isubj,:),'k')
        hold on;
        plot(xlims,simerror{ipriority}(isubj,:),'Color',aspencolors('booger'));
        defaultplot
        if ipriority == 1; title('euclidean error'); end
        
        if (expnumber == 2)
        % histogram of disc size
        datacounts = hist(data{isubj}{ipriority}(:,2),xlims);
        simdatacounts = hist(preddata{isubj}{ipriority}(:,2),xlims);
        discsize{ipriority}(isubj,:) = datacounts./sum(datacounts);
        simdiscsize{ipriority}(isubj,:)  = simdatacounts./sum(simdatacounts);
        
        subplot(3,2,2*ipriority)
        plot(xlims,discsize{ipriority}(isubj,:),'k')
        hold on;
        plot(xlims,simdiscsize{ipriority}(isubj,:),'Color',aspencolors('booger'));
        defaultplot
        if ipriority == 1; title('disc size'); end
        end
    end
    pause;
end

%% group plot
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

figure;
colorMat = {'r','b','k'};
for ipriority = 1:nPriorities
    
    % error
    if (expnumber == 2)
    subplot(3,3,ipriority)
    else
        subplot(1,3,ipriority)
    end
    fill([xlims fliplr(xlims)],[meansimerror{ipriority}-semsimerror{ipriority}...
        fliplr(meansimerror{ipriority}+semsimerror{ipriority})],colorMat{ipriority},'EdgeColor','none','FaceAlpha',0.4);
    hold on;
    errorbar(xlims,meanerror{ipriority},semerror{ipriority},'Color','k');
    defaultplot
    axis([0 10 0 0.6])
    if ipriority == 1; title('euclidean error'); end
    
    if (expnumber == 2)
        % discsize
        subplot(3,3,3+ipriority)
        fill([xlims fliplr(xlims)],[meansimdiscsize{ipriority}-semsimdiscsize{ipriority}...
            fliplr(meansimdiscsize{ipriority}+semsimdiscsize{ipriority})],colorMat{ipriority},'EdgeColor','none','FaceAlpha',0.4);
        hold on;
        errorbar(xlims,meandiscsize{ipriority},semdiscsize{ipriority},'Color','k');
        defaultplot
        axis([0 10 0 0.6])
        if ipriority == 1; title('disc size'); end
    end
end

%% quantile correlation plot per subject

nQuants = 6;
for isubj = 1:11
    figure;
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
        
        subplot(3,1,ipriority)
        plot(meanquanterror{ipriority}(isubj,:),meanquantdiscsize{ipriority}(isubj,:),'k');
        hold on;
        plot(meanquantsimerror{ipriority}(isubj,:),meanquantsimdiscsize{ipriority}(isubj,:),'Color',aspencolors('booger'));
    end
    pause;
end

%% group plot

meanmeanquanterror = cellfun(@mean,meanquanterror,'UniformOutput',false);
meanmeanquantdiscsize = cellfun(@mean,meanquantdiscsize,'UniformOutput',false);
semmeanquantdiscsize= cellfun(@(x) std(x)./sqrt(size(x,1)),meanquantdiscsize,'UniformOutput',false);
meanmeanquantsimerror = cellfun(@nanmean,meanquantsimerror,'UniformOutput',false);
meanmeanquantsimdiscsize = cellfun(@mean,meanquantsimdiscsize,'UniformOutput',false);
semmeanquantsimdiscsize= cellfun(@(x) std(x)./sqrt(size(x,1)),meanquantsimdiscsize,'UniformOutput',false);

% figure;
colorMat = {'r','b','k'};
for ipriority = 1:nPriorities
        subplot(3,3,6+ipriority)
    hold on
    plot_summaryfit(meanmeanquantsimerror{ipriority},[],[],meanmeanquantsimdiscsize{ipriority},...
    semmeanquantsimdiscsize{ipriority},[],colorMat{ipriority})
plot_summaryfit(meanmeanquanterror{ipriority},meanmeanquantdiscsize{ipriority},semmeanquantdiscsize{ipriority},...
    [],[],'k')

end

%% % % % % % % % % % % % % % % % % % % % % % % 
%           SIMULATE DATA
% % % % % % % % % % % % % % % % % % % % % % % 

clear all
expnumber = 1;
imodel = 3;
nSubj = 10;

% switch model
%     case 1 
%         logflag = logical([1 1 0]);
%     case 2
%         logflag = logical([1 1 0 0 0]);
% end

filename = ['fits/exp' num2str(expnumber) '/'];
load([filename 'fits_model' num2str(imodel) '.mat'])
% ML_parameters(logflag) = log(ML_parameters(logflag));
MU = mean(ML_parameters);
SIGMA = cov(ML_parameters);

simtheta = mvnrnd(MU,SIGMA,nSubj);
simtheta = abs(simtheta); % hacky way to enforce positive parameter values

% simtheta(:,logflag) = exp(simtheta(:,logflag));

nTrials = [250 120 70]; % mean number of trials across actual participants
for isubj = 1:10
    isubj
    simdata{isubj} = simulate_data(imodel,expnumber, simtheta(isubj,:),nTrials);
end
save([filename 'simdata_model' num2str(imodel) '.mat'],'simdata','simtheta')
