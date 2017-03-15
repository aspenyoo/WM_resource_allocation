%% understanding Jbar, tau, and sigma

Jbar = 3;
tau = 1.5;
beta = 0.01;

[JVec] = loadvar({'JVec',Jbar,tau});
% JVec = linspace(1,2,100);
Jpdf = gampdf(JVec,Jbar/tau,tau);
% Jpdf = Jpdf./qtrapz(Jpdf); % normalize

p_r = calc_pdf_r(beta,Jbar);

figure; 
subplot(2,1,1)
plot(JVec,Jpdf,'k-'); defaultplot
subplot(2,1,2)
plot(p_r,'k-'); defaultplot



%% ==========================================================
%             STUFF WITH ACTUAL DATA!
% ===========================================================
% 02.5.2017

% load in data file
load('group_data.mat')

% saving subjid, priority, radius, response_x, response_y, target_theta
usefuldata = group_data(:,[1 2 4 6 7 13]);

% deleting nans
idx = logical(sum(isnan(usefuldata),2)); % idxs of nans 
usefuldata(idx,:) = [];

% rename useful columns of group_data.mat
data_subj = usefuldata(:,1);
data_priority = usefuldata(:,2);
data_radius = usefuldata(:,3);
response_x = usefuldata(:,4);
response_y = usefuldata(:,5);
target_x = 10*cosd(usefuldata(:,6)); % x location
target_y = 10*sind(usefuldata(:,6)); % y location

% errors
error_x = response_x - target_x;
error_y = response_y - target_y;
error_distance = sqrt(error_x.^2 + error_y.^2);

% deleting trials with technical malfunctions
idx = (error_distance >= 10) | (data_radius >= 10);
error_distance(idx) = [];
data_priority(idx) = [];
data_subj(idx) = [];
data_radius(idx) = [];
error_x(idx) = [];
error_y(idx) = [];

% priority and subject numbers
priorityVec = sort(unique(data_priority),'descend'); % priority condition
nPriorities = length(priorityVec);
titleVec = {'high','med','low'};
subjVec = unique(data_subj);
nSubj = length(subjVec);

% data = cell(1,nSubj);
% for isubj = 1:nSubj
%     subjnum = subjVec(isubj);
%     
%     data{isubj} = cell(1,nPriorities);
%     for ipriority = 1:nPriorities
%         priority = priorityVec(ipriority);
%         idx = (data_subj == subjnum) & (data_priority == priority);
%         nTrials = sum(idx);
%         
%         data{isubj}{ipriority} = nan(nTrials,2);
%         data{isubj}{ipriority}(:,1) = error_distance(idx);
%         data{isubj}{ipriority}(:,2) = data_radius(idx);
%     end
% end

% save('cleandata.mat','data')

%% % % % % % % % % % % % % % % % % % % % % % % %
%     GET ML PARAMETER ESTIMATES 
% % % % % % % % % % % % % % % % % % % % % % % % 

model = 2;
nSubj = 10;
fakedata = 1;
if (fakedata)
    pretxt = 'paramrecov';
else
    pretxt = 'fits';
end

switch model
    case 1
        nParams = 3;
    case 2 
        nParams = 5;
end

bfp = nan(nSubj,nParams);
nLL = nan(1,nSubj);
for isubj = 1:nSubj
    isubj
    load([pretxt '_model' num2str(model) '_subj' num2str(isubj) '.mat'])
    blah = ML_parameters(nLLVec == min(nLLVec),:);
    bfp(isubj,:) = blah(1,:);
    nLL(isubj) = min(nLLVec);
end
ML_parameters = bfp;
nLLVec = nLL;
save([pretxt '_model' num2str(model) '.mat'],'ML_parameters','nLLVec')

%% parameter recovery plot

clear all
model = 1;

load(['paramrecov_model' num2str(model) '.mat'])
bfp = ML_parameters;

load(['simdata_model' num2str(model) '.mat'])
nParams = size(bfp,2);

for iparam = 1:nParams
    subplot(2,3,iparam);
    plot(bfp(:,iparam),simtheta(:,iparam),'ko'); hold on; 
    plot([min([bfp(:,iparam);simtheta(:,iparam)]),max([bfp(:,iparam);simtheta(:,iparam)])],...
        [min([bfp(:,iparam);simtheta(:,iparam)]),max([bfp(:,iparam);simtheta(:,iparam)])],'k-')
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
%     nLLVec3(isubj) = calc_nLL(imodel,[log(simtheta(isubj,1:2)) simtheta(isubj,3:end)],simdata{isubj});
end

[nLLVec; nLLVec2; nLLVec3; nLLVec4]

%% make nLL landscape

clear all
imodel = 1;
isubj = 6;

load(['simdata_model' num2str(imodel) '.mat'])
load(['paramrecov_model' num2str(imodel) '.mat'])

thetas = sort([simtheta(isubj,:); ML_parameters(isubj,:)]);

JbartotalVec = linspace(log(thetas(1,1).*0.5),log(thetas(2,1)*2),11);
tauVec = linspace(max(log([ 1e-5 thetas(1,2).*0.5])),max(log([1e-5 thetas(2,2)*2])),11);
% betaVec = linspace(thetas(3).*0.5,thetas(3)*1.5,11);

beta = simtheta(isubj,3);
for iJ = 1:11
    Jbartotal = JbartotalVec(iJ)
    
    for itau = 1:11
        tau = tauVec(itau);
        
        nLLMat(iJ,itau) = calc_nLL(imodel,[Jbartotal tau beta],simdata{isubj});
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
    

%% optimal pVec for model 1

load(['fits_model' num2str(model) '.mat'])
nSubj = 11;

pMat = nan(nSubj,3);
for isubj = 1:nSubj
    
    pMat(isubj,:) = calc_optimal_pVec(ML_parameters(isubj,:));
end

pMat

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

nParamVec = nan(1,2);
for imodel = 1:2
    load(['fits_model' num2str(imodel) '.mat'])
    nLL.(['model' num2str(imodel)]) = nLLVec;
    nParamVec(imodel) = size(ML_parameters,2);
    AIC.(['model' num2str(imodel)]) = 2*nLLVec + 2*nParamVec(imodel);
end

modcompidx = 1;
comparison = structfun(@(x) x-AIC.(['model' num2str(modcompidx)]),AIC,'UniformOutput',false);
comparison = comparison.model2;
meancomp = mean(comparison);
semcomp = std(comparison)/sqrt(length(comparison));

figure;
fill([0 12 12 0],[meancomp-semcomp meancomp-semcomp meancomp+semcomp meancomp+semcomp],...
    0.7*ones(1,3),'EdgeColor','none'); 
hold on;
bar(comparison,'k')

defaultplot
set(gca,'XTick',[],'XTickLabel',[])
ylabel('\Delta AIC (favoring optimal model)')


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

nPriorities = 3;
model = 1;
isubj = 1; % fitted data isubj = actual data isubj - 3
load('cleandata.mat','data')

% get ML parameter estimate for isubj
load(['fits_model' num2str(model) '.mat'])
Theta = ML_parameters(isubj,:); 

nTrials = nan(1,3);
for ipriority = 1:nPriorities
    nTrials(ipriority) = size(data{isubj}{ipriority},1);
end

simdata = simulate_data(model,Theta,nTrials);
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
    
    % plot discsize
    subplot(3,2,2*ipriority)
    hist([data{isubj}{ipriority}(:,2) simdata{ipriority}(:,2)])
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

nPriorities = 3;
model = 2;
nTrials = 1e5*ones(1,3); % how many trials to simulate per priority
nSubj = 11;
load('cleandata.mat','data')

% get ML parameter estimate for isubj
load(['fits_model' num2str(model) '.mat'])

simdata = cell(1,nSubj);
for isubj = 1:nSubj
    isubj
    Theta = ML_parameters(isubj,:);
    simdata{isubj} = simulate_data(model,Theta,nTrials);
end

%% histograms per subjects
xlims = linspace(0,10,11);

for isubj = 1:11
    figure;
    for ipriority = 1:nPriorities
        
        % histogram of euclidean error
        datacounts = hist(data{isubj}{ipriority}(:,1),xlims);
        simdatacounts = hist(simdata{isubj}{ipriority}(:,1),xlims);
        error{ipriority}(isubj,:) = datacounts./sum(datacounts);
        simerror{ipriority}(isubj,:) = simdatacounts./sum(simdatacounts);
        
        subplot(3,2,2*ipriority-1)
        plot(xlims,error{ipriority}(isubj,:),'k')
        hold on;
        plot(xlims,simerror{ipriority}(isubj,:),'Color',aspencolors('booger'));
        defaultplot
        if ipriority == 1; title('euclidean error'); end
        
        % histogram of disc size
        datacounts = hist(data{isubj}{ipriority}(:,2),xlims);
        simdatacounts = hist(simdata{isubj}{ipriority}(:,2),xlims);
        discsize{ipriority}(isubj,:) = datacounts./sum(datacounts);
        simdiscsize{ipriority}(isubj,:)  = simdatacounts./sum(simdatacounts);
        
        subplot(3,2,2*ipriority)
        plot(xlims,discsize{ipriority}(isubj,:),'k')
        hold on;
        plot(xlims,simdiscsize{ipriority}(isubj,:),'Color',aspencolors('booger'));
        defaultplot
        if ipriority == 1; title('disc size'); end
    end
    pause;
end

%% group plot
meanerror = cellfun(@mean,error,'UniformOutput',false);
semerror = cellfun(@(x) std(x)./sqrt(size(x,1)),error,'UniformOutput',false);
meandiscsize = cellfun(@mean,discsize,'UniformOutput',false);
semdiscsize = cellfun(@(x) std(x)./sqrt(size(x,1)),discsize,'UniformOutput',false);

meansimerror = cellfun(@mean,simerror,'UniformOutput',false);
semsimerror = cellfun(@(x) std(x)./sqrt(size(x,1)),simerror,'UniformOutput',false);
meansimdiscsize = cellfun(@mean,simdiscsize,'UniformOutput',false);
semsimdiscsize = cellfun(@(x) std(x)./sqrt(size(x,1)),simdiscsize,'UniformOutput',false);

figure
for ipriority = 1:nPriorities
    
    % error
    subplot(3,3,3*ipriority-2)
    fill([xlims fliplr(xlims)],[meansimerror{ipriority}-semsimerror{ipriority}...
        fliplr(meansimerror{ipriority}+semsimerror{ipriority})],aspencolors('booger'),'EdgeColor','none');
    hold on;
    errorbar(xlims,meanerror{ipriority},semerror{ipriority},'Color','k');
    defaultplot
    axis([0 10 0 0.6])
    if ipriority == 1; title('euclidean error'); end
    
    % discsize
    subplot(3,3,3*ipriority-1)
    fill([xlims fliplr(xlims)],[meansimdiscsize{ipriority}-semsimdiscsize{ipriority}...
        fliplr(meansimdiscsize{ipriority}+semsimdiscsize{ipriority})],aspencolors('booger'),'EdgeColor','none');
    hold on;
    errorbar(xlims,meandiscsize{ipriority},semdiscsize{ipriority},'Color','k');
    defaultplot
    axis([0 10 0 0.6])
    if ipriority == 1; title('disc size'); end
end

%% quantile correlation plot per subject

nQuants = 6;
for isubj = 1:11
    figure;
    for ipriority = 1:nPriorities
        currdata = data{isubj}{ipriority}(:,1);
        [currdata,idx] = sort(currdata);
        quantVec = round(linspace(0,length(currdata),nQuants+1));
        
        currsimdata = simdata{isubj}{ipriority}(:,1);
        [currsimdata,simidx] = sort(currsimdata);
        simquantVec = round(linspace(0,length(currsimdata),nQuants+1));
        for iquant = 1:nQuants
            meanquanterror{ipriority}(isubj,iquant) = mean(currdata(quantVec(iquant)+1:quantVec(iquant+1)));
            meanquantdiscsize{ipriority}(isubj,iquant) = mean(data{isubj}{ipriority}(idx(quantVec(iquant)+1:quantVec(iquant+1)),2));
            
            meanquantsimerror{ipriority}(isubj,iquant) = mean(currsimdata(simquantVec(iquant)+1:simquantVec(iquant+1)));
            meanquantsimdiscsize{ipriority}(isubj,iquant) = mean(simdata{isubj}{ipriority}(simidx(simquantVec(iquant)+1:simquantVec(iquant+1)),2));
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
meanmeanquantsimerror = cellfun(@mean,meanquantsimerror,'UniformOutput',false);
meanmeanquantsimdiscsize = cellfun(@mean,meanquantsimdiscsize,'UniformOutput',false);
semmeanquantsimdiscsize= cellfun(@(x) std(x)./sqrt(size(x,1)),meanquantsimdiscsize,'UniformOutput',false);

% figure;
colorMat = {'k','r','b'};
for ipriority = 1:nPriorities
        subplot(3,3,3*ipriority)
    hold on
    plot_summaryfit(meanmeanquantsimerror{ipriority},[],[],meanmeanquantsimdiscsize{ipriority},...
    semmeanquantsimdiscsize{ipriority},[],aspencolors('booger'))
plot_summaryfit(meanmeanquanterror{ipriority},meanmeanquantdiscsize{ipriority},semmeanquantdiscsize{ipriority},...
    [],[],'k')

end

%% % % % % % % % % % % % % % % % % % % % % % % 
%           SIMULATE DATA
% % % % % % % % % % % % % % % % % % % % % % % 

clear all
model = 2;
nSubj = 10;

% switch model
%     case 1 
%         logflag = logical([1 1 0]);
%     case 2
%         logflag = logical([1 1 0 0 0]);
% end

load(['fits_model' num2str(model) '.mat'])
% ML_parameters(logflag) = log(ML_parameters(logflag));
MU = mean(ML_parameters);
SIGMA = cov(ML_parameters);

simtheta = mvnrnd(MU,SIGMA,nSubj);
simtheta = abs(simtheta); % hacky way to parameters are positive

% simtheta(:,logflag) = exp(simtheta(:,logflag));

nTrials = [250 120 70]; % mean number of trials across actual participants
for isubj = 1:nSubj
    isubj
    simdata{isubj} = simulate_data(model,simtheta(isubj,:),nTrials);
end
save(['simdata_model' num2str(model) '.mat'],'simdata','simtheta')