
%% ====================================
%              MAIN FIGURES
%  ====================================
%
% fig 1b: exp 1 main behavioral results
% fig 2b: exp 1 modeling results
% fig 2c/4b: exp 1/2 ternary plot
% fig 2d/4c: exp 1/2 model comparison
% fig 3b: exp 2 main behavioral results
% fig 4a: exp 2 modeling results

%% FIG 1B: EXP 1 MAIN BEHAVIORAL RESULTS

clear all; close all

load('exp1_cleandata.mat','data')
nSubj = length(data);

priorityVec = [0.1 0.3 0.6];
nPriorities = 3;
colorMat = [0 0 0; 0 0 1; 1 0 0];

figure('pos',[500 500 300 400])
hold on

% get error data
dataerrorVec = nan(nSubj,nPriorities);
for ipriority = 1:nPriorities
    for isubj = 1:nSubj
        dataerrorVec(isubj,ipriority) = mean(data{isubj}{ipriority}(:,1));
    end
end
% M and SEM
dataMerror = mean(dataerrorVec)
dataSEMerror = std(dataerrorVec)./sqrt(nSubj)

% plot
dx = 0.2;
defaultplot;
axis([0.5 3.5 0 2.2])
set(gca,'YTick',0:.2:2.2);
for ipriority = 1:nPriorities
    errorb(ipriority,dataMerror(nPriorities+1-ipriority),dataSEMerror(nPriorities+1-ipriority),...
        'LineWidth',1,'Color',colorMat(ipriority,:));
end
xlabel('probe probability'); ylabel('error (dva)');
set(gca,'XTick',1:3,'XTickLabel',priorityVec,...
    'YTick',0:0.2:2.2,'YTickLabel',{0,'','','','',1,'','','','',2,''})


%% plot primary and final saccade


M_pri = fliplr([1.8482 2.1755 2.5740]);
sem_pri = fliplr([0.0942 0.1409 0.1930]);

M_final = fliplr([1.2835    1.5830    1.8749]);
sem_final = fliplr([0.0703    0.1170    0.1711]);

% is this as a function of primary or final saccade?
M_rt = fliplr([463.7298  469.4223  478.5309]);
sem_rt = fliplr([14.6151  16.6736  17.6714]);

priorityVec = fliplr(priorityVec);

figure;
subplot(1,3,1)
defaultplot;
axis([0.5 3.5 0 3])
errorb(M_pri,sem_pri);
xlabel('probe probability'); ylabel('error (dva)');
set(gca,'XTick',1:3,'XTickLabel',priorityVec,...
    'YTick',0:0.5:3)

subplot(1,3,2)
defaultplot;
axis([0.5 3.5 0 3])
errorb(M_final,sem_final);
xlabel('probe probability'); ylabel('error (dva)');
set(gca,'XTick',1:3,'XTickLabel',priorityVec,...
    'YTick',0:0.5:3)

subplot(1,3,3)
defaultplot;
axis([0.5 3.5 0 600])
errorb(M_rt,sem_rt);
xlabel('probe probability'); ylabel('RT (msec)');
set(gca,'XTick',1:3,'XTickLabel',priorityVec)


%% FIG 2B: EXP 1 MODELING RESULTS

clear all; close all

fixedrisk = [];
colorMat = [1 0 0; 0 0 1; 0 0 0; 0 1 0];
modelorderVec = [3 2 4]; % the order the models should be plotted
nModels = length(modelorderVec); % how many models you want to plot now
modelnameVec = {'max points','flexible','proportional','min error'};

nSubj = 14;
nPriorities = 3;
xlims = linspace(0,10,16); % for histograms
xdiff = diff(xlims(1:2));
filepath = 'fits/priority/exp1/';

figure;
ha = tight_subplot(1,nModels,[.03 .03],[.2 .07],[.06 .01]);
set(gcf,'Position',[28 504 800 236])

% load and plot predicted data for each model
for imodel = 1:nModels
    model = modelorderVec(imodel);
    
    % load data
    load([filepath 'modelpred_exp1_model' num2str(model) fixedrisk '.mat'],'preddata')
    
    % get model predictions
    for isubj = 1:nSubj
        for ipriority = 1:nPriorities
            % histogram of euclidean error
            simdatacounts = hist(preddata{isubj}{ipriority}(:,1),xlims);
            simerror{ipriority}(isubj,:) = simdatacounts./sum(simdatacounts);
        end
    end
    
    meansimerror = cellfun(@mean,simerror,'UniformOutput',false);
    semsimerror = cellfun(@(x) std(x)./sqrt(size(x,1)),simerror,'UniformOutput',false);
    
    % plot error (and circle size) model predictions
    for ipriority = 1:nPriorities
        
        axes(ha(imodel))
        fill( [0 xlims+xdiff fliplr(xlims)+xdiff 0],[0 meansimerror{ipriority}-semsimerror{ipriority}...
            fliplr(meansimerror{ipriority}+semsimerror{ipriority}) 0],colorMat(ipriority,:),'EdgeColor','none','FaceAlpha',0.4);
        hold on;
        title(modelnameVec{model})
    end
end

% ============= PLOT REAL DATA ======================

% load data
load('exp1_cleandata.mat','data')

% histograms per subjects
for isubj = 1:nSubj
    for ipriority = 1:nPriorities
        
        % histogram of euclidean error
        datacounts = hist(data{isubj}{ipriority}(:,1),xlims);
        error{ipriority}(isubj,:) = datacounts./sum(datacounts);
    end
end

% =========== group plot =====================

for imodel = 1:nModels
    
    meanerror = cellfun(@mean,error,'UniformOutput',false);
    semerror = cellfun(@(x) std(x)./sqrt(size(x,1)),error,'UniformOutput',false);
    
    for ipriority = 1:nPriorities
        
        
        % =========== ERROR ============
        
        axes(ha(imodel))
        hold on;
        errorbar(xlims+xdiff,meanerror{ipriority},semerror{ipriority},'Color',colorMat(ipriority,:),'LineStyle', 'none','LineWidth',1);
        defaultplot
        axis([0 10 0 0.4])
        
        xlabel('error (dva)','FontSize',16)
        set(ha(imodel),'YTick',0:0.1:0.4,'XTick',0:10);
        
        if imodel == 1
            ylabel('proportion','FontSize',16)
        else
            set(ha(imodel),'YTickLabel','');
        end
    end
end


%% FIG 2C/4B: EXP 1/2 TERNARY PLOT

% clear all; close all
expnumber = 1; % experiment number. change to toggle between experiments

% load data
load(sprintf('fits/priority/exp%d/fits_model2.mat',expnumber))
pMat = ML_parameters(:,end-1:end);
pMat(:,3) = 1-sum(pMat,2);

% axis
figure;
[h,hg,htick]=terplot;
c1 = [1 0.5 1/3];
c2 = [0 0.5 1/3];
hlabels=terlabel('high','medium','low');
set(h,'LineWidth',1)
colorMat = [1 0 0; 0 0 1; 0 0 0];

% define colors
axis1 = colorMat(2,:);
axis2 = colorMat(1,:);
axis3 = colorMat(3,:);
grey = 0.7*ones(1,3);

% make patch showing monotonic area
x=0.5-c1*cos(pi/3)+c2/2;
y=0.866-c1*sin(pi/3)-c2*cot(pi/6)/2;
patch('Faces',[1 2 3],'Vertices',[x' y'],'FaceColor',grey,'FaceAlpha',0.3,'EdgeColor','none');

% change the color of the grid lines
set(hg(:,1),'color',axis1)
set(hg(:,2),'color',axis2)
set(hg(:,3),'color',axis3)

% make 0.6 0.3 0.1 lines noticeable
set(hg(3,1),'LineStyle','-')
set(hg(6,2),'LineStyle','-')
set(hg(1,3),'LineStyle','-')

% modify the label size and color
set(hlabels,'fontsize',12)
set(hlabels(1),'color',axis1)
set(hlabels(2),'color',axis2)
set(hlabels(3),'color',axis3)

% modify the tick label colors
set(htick(:,1),'color',axis1,'linewidth',3)
set(htick(:,2),'color',axis2,'linewidth',3)
set(htick(:,3),'color',axis3,'linewidth',3)

% plot data
indigo = [90 90 190]./255;
hter=ternaryc(pMat(:,1),pMat(:,2),pMat(:,3));
set(hter,'marker','o','markerfacecolor',indigo,'markersize',8,'markeredgecolor',indigo)

% plot me points
redish = [221 116 88]./255;
hter=ternaryc(pMat_me(:,1),pMat_me(:,2),pMat_me(:,3));
set(hter,'marker','o','markerfacecolor',redish,'markersize',8,'markeredgecolor',redish)

% plot connecting lines
c1 = pMat(:,1);
c2 = pMat(:,2);
x1=0.5-c1*cos(pi/3)+c2/2;
y1=0.866-c1*sin(pi/3)-c2*cot(pi/6)/2;
c1 = pMat_me(:,1);
c2 = pMat_me(:,2);
x2=0.5-c1*cos(pi/3)+c2/2;
y2=0.866-c1*sin(pi/3)-c2*cot(pi/6)/2;
hold on;
plot([x1 x2]',[y1 y2]','--k')

%% FIG 2D/4C: EXP 1/2 MODEL COMPARISON

clear all
expnumber = 2; % experiment number. change to toggle between experiments

switch expnumber
    case 1
        modVec = [2 3 4]; % [flex prop min_error]
    case 2
        modVec = [2 3 4 1]; % [flex prop min_error max_points]
end
nModels = length(modVec);
modcompidx = 2; % relative to flexible model

filename = ['exp' num2str(expnumber) '_cleandata.mat'];
load(filename)

nTrials = sum(nTrials,2);
nSubj = length(data);

filepath = ['fits/priority/exp' num2str(expnumber) '/'];

nLLMat = nan(nModels,nSubj);
for imodel = 1:nModels
    imodel
    modidx = modVec(imodel);
    
    load([filepath 'fits_model' num2str(modidx) '.mat'])
    nLLMat(imodel,:) = nLLVec;
    nParamVec(imodel) = size(ML_parameters,2);
end

% get AICc and BIC
modelcomparisons = cell(1,2);
[modelcomparisons{2}, modelcomparisons{1}]= modcomp(nLLMat',nParamVec,nTrials);

% labels and index stuff
modlabels = {'max points','flexible','prop','min error'};
modcomplabel = modlabels{modcompidx};
modlabels = modlabels(modVec(modVec ~= modcompidx));

nSubj = length(nLLVec);
subtractything = 2-expnumber;

fillx = 0.4;
titleVec = {'\Delta AICc','\Delta BIC'};

figure;
for iplot = 1:2;
    subplot(1,2,iplot)
    comparison = bsxfun(@minus,modelcomparisons{iplot},modelcomparisons{iplot}(:,modVec == modcompidx));
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
    
    xthing = [];
    if nModels == 2
        fill([0 nSubj+1 nSubj+1 0],medCI([1 1 2 2]),...
            0.8*ones(1,3),'EdgeColor','none'); hold on;
        plot([0 nSubj+1],[mediancomp mediancomp],'Color',[0.1 0.1 0.1])
        set(gca,'XTick',[],'XTickLabel',modlabels)
    else
        for imodel = 1:(nModels-1)
            fill([imodel-fillx imodel+fillx imodel+fillx imodel-fillx],medCI([imodel imodel imodel+nModels-1 imodel+nModels-1]),...
                0.7*ones(1,3),'EdgeColor','none'); hold on;
            plot([imodel-fillx imodel+fillx],mediancomp(imodel)*ones(1,2),'Color',[0.1 0.1 0.1])
            set(gca,'XTick',[],'XTick',1:3,'XTickLabel',modlabels)
            xthing = [xthing imodel*ones(nSubj,1)];
        end
    end
    
    for isubj = 1:nSubj
        plot(xthing(isubj,:),comparison(isubj,:),'k.','MarkerSize',24)
    end
    plot([0.5 nModels-0.5],[0 0],'k-') % 0 axis line
    
    title(titleVec{iplot})
    defaultplot
    
end

%% FIG 3B: EXP 2 MAIN BEHAVIORAL RESULTS

clear all

exppriorityVec = [0.6 0.3 0.1];
nPriorities = length(exppriorityVec);
colorMat = [0 0 0; 0 0 1; 1 0 0];

expnumber = 2;
load('exp2_cleandata.mat','data')
nSubj = 11;

[dataerrorVec, datadiscsizeVec] = deal(nan(nSubj,nPriorities));
for ipriority = 1:nPriorities
    
    for isubj = 1:nSubj
        dataerrorVec(isubj,ipriority) = mean(data{isubj}{ipriority}(:,1));
        datadiscsizeVec(isubj,ipriority) = mean(data{isubj}{ipriority}(:,2));
    end
    
end

%  ========= calculate M, SEM across subjects ========
dataMerror = mean(dataerrorVec);
dataSEMerror = std(dataerrorVec)./sqrt(nSubj);

dataMdiscsize = mean(datadiscsizeVec);
dataSEMdiscsize = std(datadiscsizeVec)/sqrt(nSubj);


% ===== plot main effect: error =======
figure('pos',[500 500 1000 300])
dx = 0.2;
subplot(1,3,1);hold on;
axis([0.5 3.5 0 2.4])
xlabel('probe probability'); ylabel('error (dva)');
set(gca,'YTick',0:0.4:2.4, 'YTickLabel', {0,'','',1.2,'','' 2.4},...
    'XTick',1:3,'XTickLabel',{'0.1','0.3','0.6'})
for ipriority = 1:nPriorities
    errorb(ipriority+dx,dataMerror(nPriorities+1-ipriority),dataSEMerror(nPriorities+1-ipriority),...
        'Color',colorMat(ipriority,:),'LineWidth',1); 
end
defaultplot;

% ===== plot main effect: radius of circle wager ======
subplot(1,3,2); hold on;
axis([0.5 3.5 0 4])
for ipriority = 1:nPriorities
    errorb(ipriority+dx,dataMdiscsize(nPriorities+1-ipriority),dataSEMdiscsize(nPriorities+1-ipriority),...
        'Color',colorMat(ipriority,:),'LineWidth',1); 
end
xlabel('probe probability'); ylabel('circle radius, r')
set(gca,'YTick',0:4,'YTickLabel',{0 '', 2,'', 4},...
    'XTick',1:3,'XTickLabel',{'0.1','0.3','0.6'})
defaultplot;

% ========== plot correlation between error and circle radius ===========

% group data into quantiles
nQuants = 6;
[meanquanterror, meanquantdiscsize] = deal(cell(1,nPriorities));
for isubj = 1:11
    for ipriority = 1:nPriorities
        currdata = data{isubj}{ipriority}(:,1);
        [currdata,idx] = sort(currdata);
        quantVec = round(linspace(0,length(currdata),nQuants+1));
        for iquant = 1:nQuants
            meanquanterror{ipriority}(isubj,iquant) = mean(currdata(quantVec(iquant)+1:quantVec(iquant+1)));
            meanquantdiscsize{ipriority}(isubj,iquant) = mean(data{isubj}{ipriority}(idx(quantVec(iquant)+1:quantVec(iquant+1)),2));
        end
    end
end

% plot means \pm SEMs of quantiles
meanmeanquanterror = cellfun(@mean,meanquanterror,'UniformOutput',false);
meanmeanquantdiscsize = cellfun(@mean,meanquantdiscsize,'UniformOutput',false);
semmeanquantdiscsize= cellfun(@(x) std(x)./sqrt(size(x,1)),meanquantdiscsize,'UniformOutput',false);

subplot(1,3,3); hold on
axis([0 6 0 4])
for ipriority = 1:nPriorities
    errorbar(meanmeanquanterror{nPriorities+1-ipriority},meanmeanquantdiscsize{nPriorities+1-ipriority},...
        semmeanquantdiscsize{nPriorities+1-ipriority},'Color',colorMat(ipriority,:),'LineWidth',1); 
    errorb(meanmeanquanterror{nPriorities+1-ipriority},meanmeanquantdiscsize{nPriorities+1-ipriority},...
        semmeanquantdiscsize{nPriorities+1-ipriority},'Color',colorMat(ipriority,:),'LineWidth',1); 
end
ylabel('circle radius, r');
set(gca,'XTick',[0 3 6],'XTickLabel',[0,3,6],...
    'YTick',0:4,'YTickLabel',{0,'',2,'',4});
xlabel('error (dva)');
defaultplot


%% FIG 4A: EXP 2 MODELING RESULTS

clear all; close all
colorMat = [1 0 0; 0 0 1; 0 0 0; 0 1 0];
modelorderVec = [3 2 4 1]; % the order the models should be plotted
nModels = length(modelorderVec); % how many models you want to plot now
modelnameVec = {'Max Points','Flexible','Proportional','Min Error'};

nSubj = 11;

nPriorities = 3;
xlims = linspace(0,10,16); % for histograms
xdiff = diff(xlims(1:2));

figure;
% top down, left right, bottom and top margins, left and right margins
ha = tight_subplot(3,nModels,{[.1 .1],[.03 .03 .03]},[.08 .03],[.05 .03]);
set(gcf,'Position',[28 504 800 500])


% =========== PLOT MODEL PREDICTIONS =============

for imodel = 1:nModels
    model = modelorderVec(imodel);
    
    % load data
    load(sprintf('fits/priority/exp2/modelpred_exp2_model%d.mat',model),'preddata')
    
    % get model predictions
    for isubj = 1:nSubj
        for ipriority = 1:nPriorities
            % histogram of euclidean error
            simdatacounts = hist(preddata{isubj}{ipriority}(:,1),xlims);
            simerror{ipriority}(isubj,:) = simdatacounts./sum(simdatacounts);
            
            % histogram of circle size
            simdatacounts = hist(preddata{isubj}{ipriority}(:,2),xlims);
            simcirclesize{ipriority}(isubj,:)  = simdatacounts./sum(simdatacounts);
        end
    end
    
    meansimerror = cellfun(@mean,simerror,'UniformOutput',false);
    semsimerror = cellfun(@(x) std(x)./sqrt(size(x,1)),simerror,'UniformOutput',false);
    
    
    meansimcirclesize = cellfun(@mean,simcirclesize,'UniformOutput',false);
    semsimcirclesize = cellfun(@(x) std(x)./sqrt(size(x,1)),simcirclesize,'UniformOutput',false);
    
    
    % plot error (and circle size) model predictions
    for ipriority = 1:nPriorities
        
        % error
        axes(ha(imodel))
        fill( [0 xlims+xdiff fliplr(xlims)+xdiff 0],[0 meansimerror{ipriority}-semsimerror{ipriority}...
            fliplr(meansimerror{ipriority}+semsimerror{ipriority}) 0],colorMat(ipriority,:),'EdgeColor','none','FaceAlpha',0.4);
        hold on;
        title(modelnameVec{model})
        
        % circle size
        axes(ha(nModels+imodel))
        fill([0 xlims fliplr(xlims) 0]+xdiff,[0 meansimcirclesize{ipriority}-semsimcirclesize{ipriority}...
            fliplr(meansimcirclesize{ipriority}+semsimcirclesize{ipriority}) 0],colorMat(ipriority,:),'EdgeColor','none','FaceAlpha',0.4);
        hold on;
    end
    
    % get quantile information per subject
    nQuants = 6;
    for isubj = 1:nSubj
        
        for ipriority = 1:nPriorities
            currsimdata = preddata{isubj}{ipriority}(:,1);
            [currsimdata,simidx] = sort(currsimdata);
            simquantVec = round(linspace(0,length(currsimdata),nQuants+1));
            
            for iquant = 1:nQuants
                meanquantsimerror{ipriority}(isubj,iquant) = mean(currsimdata(simquantVec(iquant)+1:simquantVec(iquant+1)));
                meanquantsimdiscsize{ipriority}(isubj,iquant) = mean(preddata{isubj}{ipriority}(simidx(simquantVec(iquant)+1:simquantVec(iquant+1)),2));
            end
        end
    end
    
    % plot correlation     
    meanmeanquantsimerror = cellfun(@nanmean,meanquantsimerror,'UniformOutput',false);
    meanmeanquantsimdiscsize = cellfun(@mean,meanquantsimdiscsize,'UniformOutput',false);
    semmeanquantsimdiscsize= cellfun(@(x) std(x)./sqrt(size(x,1)),meanquantsimdiscsize,'UniformOutput',false);
    
    axes(ha(2*nModels+imodel)); hold on
    for ipriority = 1:nPriorities
        x = meanmeanquantsimerror{ipriority};
        meanRadius = meanmeanquantsimdiscsize{ipriority};
        semRadius = semmeanquantsimdiscsize{ipriority};
        
        fill( [x fliplr(x)],[meanRadius-semRadius fliplr(meanRadius+semRadius)],...
            colorMat(ipriority,:),'EdgeColor','none','FaceAlpha',0.4);
    end
    
end

% ================= PLOT REAL DATA ======================

% load data
load('exp2_cleandata.mat','data')

% histograms per subjects
for isubj = 1:nSubj
    for ipriority = 1:nPriorities
        
        % histogram of euclidean error
        datacounts = hist(data{isubj}{ipriority}(:,1),xlims);
        error{ipriority}(isubj,:) = datacounts./sum(datacounts);
        
        % histogram of disc size
        datacounts = hist(data{isubj}{ipriority}(:,2),xlims);
        discsize{ipriority}(isubj,:) = datacounts./sum(datacounts);
    end
    
end

for imodel = 1:nModels
    
    meanerror = cellfun(@mean,error,'UniformOutput',false);
    semerror = cellfun(@(x) std(x)./sqrt(size(x,1)),error,'UniformOutput',false);
    
    meandiscsize = cellfun(@mean,discsize,'UniformOutput',false);
    semdiscsize = cellfun(@(x) std(x)./sqrt(size(x,1)),discsize,'UniformOutput',false);
    
    for ipriority = 1:nPriorities
        
        
        % =========== error ============
        axes(ha(imodel))
        hold on;
        errorbar(xlims+xdiff,meanerror{ipriority},semerror{ipriority},'Color',colorMat(ipriority,:),'LineStyle', 'none','LineWidth',1);
        defaultplot
        axis([0 10 0 0.4])
        
        xlabel('error (dva)');
        set(ha(imodel),'XTick',0:10,'YTick',0:0.1:0.4,...
            'YTickLabel','','FontSize',10);
        if (imodel == 1)
            ylabel('proportion','FontSize',14);
            set(ha(imodel),'YTickLabel',0:0.1:0.4)
        end
        
        % ========= circle size ===========
        axes(ha(nModels+imodel))
        hold on;
        errorbar(xlims+xdiff,meandiscsize{ipriority},semdiscsize{ipriority},'Color',colorMat(ipriority,:),'LineStyle','none','LineWidth',1);
        defaultplot
        axis([0 10 0 0.4])
        set(ha(nModels+imodel),'XTick',0:10,'YTick',0:0.1:0.4,...
            'FontSize',10,'YTickLabel','');
        xlabel('circle size, r','FontSize',10);
        if (imodel == 1);
            ylabel('proportion','FontSize',14);
            set(ha(nModels+imodel),'YTickLabel',0:0.1:0.4);
        end
    end
    
    % quantile binning data
    nQuants = 6;
    for isubj = 1:nSubj
        
        for ipriority = 1:nPriorities
            currdata = data{isubj}{ipriority}(:,1);
            [currdata,idx] = sort(currdata);
            quantVec = round(linspace(0,length(currdata),nQuants+1));
            
            for iquant = 1:nQuants
                meanquanterror{ipriority}(isubj,iquant) = mean(currdata(quantVec(iquant)+1:quantVec(iquant+1)));
                meanquantdiscsize{ipriority}(isubj,iquant) = mean(data{isubj}{ipriority}(idx(quantVec(iquant)+1:quantVec(iquant+1)),2));
            end
        end
    end
    
    % ================ correlation plot ====================
    meanmeanquanterror = cellfun(@mean,meanquanterror,'UniformOutput',false);
    meanmeanquantdiscsize = cellfun(@mean,meanquantdiscsize,'UniformOutput',false);
    semmeanquantdiscsize= cellfun(@(x) std(x)./sqrt(size(x,1)),meanquantdiscsize,'UniformOutput',false);
    
    axes(ha(2*nModels+imodel))
    axis([0 6 0 6])
    for ipriority = 1:nPriorities
        hold on
        errorb(meanmeanquanterror{ipriority},meanmeanquantdiscsize{ipriority},semmeanquantdiscsize{ipriority},'Color',colorMat(ipriority,:));
    end
    
    if (imodel == 1);
        ylabel('circle size, r','FontSize',14);
    end
    set(ha(2*nModels+imodel),'XTick',0:6,'YTick',0:6,'XTickLabel',0:6,'YTickLabel',0:6);
    
    xlabel('error (dva)');
end


%% ====================================
%        SUPPLEMENTARY FIGURES
%  ====================================
%
% fig 1: how optimal allocation changes with Jbar for each model
% fig 2: how optimal allocation changes with experimental probe probability

%% FIG 1:  how optimal allocation changes with Jbar for each model
% note that the colors come out wrong in this plot, but the locations are
% correct

clear all, close all
exppriorityVec = [0.6 0.3 0.1];

% ====== SET UP TERNARY PLOT ======= 
colorMat = [1 0 0; 0 0 1; 0 0 0];

% axis
figure;
[h,hg,htick]=terplot;
c1 = [1 0.5 1/3];
c2 = [0 0.5 1/3];
c3 = [0 0 1/3];
hlabels=terlabel('high','medium','low');
set(h,'LineWidth',1)

% define colors
axis1 = colorMat(2,:);
axis2 = colorMat(1,:);
axis3 = colorMat(3,:);
grey = 0.7*ones(1,3);

% change the color of the grid lines
set(hg(:,1),'color',axis1)
set(hg(:,2),'color',axis2)
set(hg(:,3),'color',axis3)

% modify the label size and color
set(hlabels,'fontsize',12)
set(hlabels(1),'color',axis1)
set(hlabels(2),'color',axis2)
set(hlabels(3),'color',axis3)

% modify the tick label colors
set(htick(:,1),'color',axis1,'linewidth',3)
set(htick(:,2),'color',axis2,'linewidth',3)
set(htick(:,3),'color',axis3,'linewidth',3)

% plot proportional model
h = ternaryc(0.6,0.3,0.1);
propcolor = [0.1 0.7 0.3];
set(h,'marker','.','markerfacecolor',propcolor,'markersize',24,'markeredgecolor',propcolor)
hold on

% CALCULATE AND PLOT MIN-ERROR AND MAX-PTS MODEL PREDICTIONS
tau = 0.4;
alpha = 1;
beta = 1;
gamma = 1;

nSamps = 10;
JbarVec = linspace(3,25,nSamps); % multiplier to use for Jbar

[pVec_MP, pVec_ME] = deal(nan(nSamps,3));
for isamp2 = 1:nSamps
    Jbar = JbarVec(isamp2)*tau + 0.01
    
    MPtheta = [Jbar tau alpha beta];
    MEtheta = [MPtheta gamma];
    
    % optimal allocation for maximizing points model
    pVec_MP(isamp2,:) = calc_pVec_maxpoints(MPtheta,exppriorityVec);
    
    % optimal allocation for minimizing error model
    pVec_ME(isamp2,:) = calc_pVec_minerror(MEtheta,exppriorityVec);
end

sz = 40;
c = linspace(1,nSamps,nSamps);

% plot max points model
x=0.5-pVec_MP(:,1).*cos(pi/3)+pVec_MP(:,2)./2;
y=0.866-pVec_MP(:,1).*sin(pi/3)-pVec_MP(:,2).*cot(pi/6)/2;
plot(x,y,'k-')
colormap('autumn')
s = scatter(x,y,sz,c,'filled')

% plot min error model
x=0.5-pVec_ME(:,1).*cos(pi/3)+pVec_ME(:,2)./2;
y=0.866-pVec_ME(:,1).*sin(pi/3)-pVec_ME(:,2).*cot(pi/6)/2;
plot(x,y,'k-')
s2 = scatter(x,y,sz,c,'filled')


%% FIG 2: how optimal allocation changes with experimental probe probability
% for Minimizing Error model for each subject in Exp 2

clear all, close all, clc

load('fits/priority/exp2/fits_model2.mat')
nSubj = size(ML_parameters,1);

% get probe probabilities that will be investigated
x = 0:0.1:1;
[X,Y] = meshgrid(x,x);
X = X(:);
Y = Y(:);
idx = (X+Y) > 1;
X(idx) = [];
Y(idx) = [];
Z = 1-X-Y;

% get x- and y-coordinates of actual values (in ternary space)
x=0.5-X.*cos(pi/3)+Y./2;
y=0.866-X.*sin(pi/3)-Y.*cot(pi/6)/2;
nSamps = length(X);

figure;
for isubj = 1:nSubj;
    isubj
    
    % parameters
    theta = ML_parameters(isubj,:);
    
    % ------ SET UP TERNARY PLOT -----
    colorMat = [1 0 0; 0 0 1; 0 0 0];
    
    % axis
    subplot(3,4,isubj);
    [h,hg,htick]=terplot;
    c1 = [1 0.5 1/3];
    c2 = [0 0.5 1/3];
    c3 = [0 0 1/3];
    hlabels=terlabel('high','medium','low');
    set(h,'LineWidth',1)
    
    % define colors
    axis1 = colorMat(2,:);
    axis2 = colorMat(1,:);
    axis3 = colorMat(3,:);
    grey = 0.7*ones(1,3);
    
    % change the color of the grid lines
    set(hg(:,1),'color',axis1)
    set(hg(:,2),'color',axis2)
    set(hg(:,3),'color',axis3)
    
    % modify the label size and color
    set(hlabels,'fontsize',12)
    set(hlabels(1),'color',axis1)
    set(hlabels(2),'color',axis2)
    set(hlabels(3),'color',axis3)
    
    % modify the tick label colors
    set(htick(:,1),'color',axis1,'linewidth',3)
    set(htick(:,2),'color',axis2,'linewidth',3)
    set(htick(:,3),'color',axis3,'linewidth',3)
    
    
    hold on
    % get optimal resource allocation for each probe probability combination
    pVec_ME = nan(nSamps,3);
    for isamp = 1:nSamps;
        exppriorityVec = [X(isamp) Y(isamp) Z(isamp)];
        pVec_ME(isamp,:) = calc_pVec_minerror(theta,exppriorityVec);
    end
    
    % get x- and y-coordinates in ternary space
    u = 0.5-pVec_ME(:,1).*cos(pi/3)+pVec_ME(:,2)./2;
    v = 0.866-pVec_ME(:,1).*sin(pi/3)-pVec_ME(:,2).*cot(pi/6)/2;
    
    % get amount of difference between optimal allocation and probe prob.
    u = u-x; 
    v = v-y;
    
    % plot
    quiver(x,y,u,v,Scale(0),'k')
    
end


%% how optimal allocation changes with experimental probe probability
% for arbitrary parameter combination for Minimizing Error model

clear all
close all
clc

% theta
theta = [5 1.6 1 1 1];

% ------ SET UP TERNARY PLOT -----
colorMat = [1 0 0; 0 0 1; 0 0 0];

% axis
figure;
[h,hg,htick]=terplot;
c1 = [1 0.5 1/3];
c2 = [0 0.5 1/3];
c3 = [0 0 1/3];
hlabels=terlabel('high','medium','low');
set(h,'LineWidth',1)

% define colors
axis1 = colorMat(2,:);
axis2 = colorMat(1,:);
axis3 = colorMat(3,:);
grey = 0.7*ones(1,3);

% change the color of the grid lines
set(hg(:,1),'color',axis1)
set(hg(:,2),'color',axis2)
set(hg(:,3),'color',axis3)

% modify the label size and color
set(hlabels,'fontsize',12)
set(hlabels(1),'color',axis1)
set(hlabels(2),'color',axis2)
set(hlabels(3),'color',axis3)

% modify the tick label colors
set(htick(:,1),'color',axis1,'linewidth',3)
set(htick(:,2),'color',axis2,'linewidth',3)
set(htick(:,3),'color',axis3,'linewidth',3)

% get probe probabilities that will be investigated
x = 0:0.1:1;
[X,Y] = meshgrid(x,x);
X = X(:);
Y = Y(:);
idx = (X+Y) > 1;
X(idx) = [];
Y(idx) = [];
Z = 1-X-Y;

hold on
% get x and y of actual values
x=0.5-X.*cos(pi/3)+Y./2;
y=0.866-X.*sin(pi/3)-Y.*cot(pi/6)/2;

nSamps = length(X);
pVec_ME = nan(nSamps,3);
for isamp = 1:nSamps;
    exppriorityVec = [X(isamp) Y(isamp) Z(isamp)];
    pVec_ME(isamp,:) = calc_pVec_minerror(theta,exppriorityVec);
end

u = 0.5-pVec_ME(:,1).*cos(pi/3)+pVec_ME(:,2)./2;
v = 0.866-pVec_ME(:,1).*sin(pi/3)-pVec_ME(:,2).*cot(pi/6)/2;

u = u-x;
v = v-y;
quiver(x,y,u,v,Scale(0),'k')

