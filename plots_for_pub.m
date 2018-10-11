
%% ====================================
%            MAIN FIGURES
%  ====================================
%
% fig 1b: exp 1 main behavioral results
% fig 2b: exp 1 modeling results
% fig 2c: exp 1 ternary plot
% fig 2d: exp 1 model comparison
% fig 3b: exp 2 main behavioral results
% fig 4a: exp 2 modeling results
% fig 4b: exp 2 ternary plot
% fig 4c: exp 2 model comparison

%% FIG 1B: EXP 1 MAIN BEHAVIORAL RESULTS

clear all; close all

load('exp1_cleandata.mat','data')
filepath = 'fits/exp1/';
nSubj = length(data);

priorityVec = [0.1 0.3 0.6];
nPriorities = 3;
colorMat = [0 0 0; 0 0 1; 1 0 0];
markerMat = {'.','x','s'};

modelnameVec = {'optimal','free','fixed'}';

figure; hold on

% get error data
dataerrorVec = nan(nSubj,nPriorities);
for ipriority = 1:nPriorities
    for isubj = 1:nSubj
        dataerrorVec(isubj,ipriority) = mean(data{isubj}{ipriority}(:,1));
    end
end
% M and SEM
dataMerror = mean(dataerrorVec);
dataSEMerror = std(dataerrorVec)./sqrt(nSubj);

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
set(gca,'XTick',1:3,'XTickLabel',{'0.1','0.3','0.6'})
title('error')

%% FIG 2B: EXP 1 MODELING RESULTS

clear all; close all

expnumber = 1;
fixedrisk = [];
colorMat = [1 0 0; 0 0 1; 0 0 0; 0 1 0];
nModelsPossible = 4; % how many total models there are
modelorderVec = [3 2 4]; % the order the models should be plotted
nModels = length(modelorderVec); % how many models you want to plot now
modelnameVec = {'max points','flexible','proportional','min error'};

nSubj = 14;
nPriorities = 3;
xlims = linspace(0,10,16); % for histograms
xdiff = diff(xlims(1:2));
filepath = 'fits/zuzanna/exp1/';

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
%% fig 2c: exp 1 ternary plot


clear all;
colorMat = [1 0 0; 0 0 1; 0 0 0];

% load data
load('fits/priority/exp1/fits_model2.mat')
nSubj = size(ML_parameters,1);
pMat = ML_parameters(:,end-1:end);
pMat(:,3) = 1-sum(pMat,2);

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

% make patch showing monotonic area
x=0.5-c1*cos(pi/3)+c2/2;
y=0.866-c1*sin(pi/3)-c2*cot(pi/6)/2;
patch('Faces',[1 2 3],'Vertices',[x' y'],'FaceColor',grey,'FaceAlpha',0.3,'EdgeColor','none');

% change the color of the grid lines
set(hg(:,1),'color',axis1)
set(hg(:,2),'color',axis2)
set(hg(:,3),'color',axis3)

% make 0.6 0.3 0.1 lines noticeable
set(hg(3,1),'LineStyle','--')
set(hg(6,2),'LineStyle','--')
set(hg(1,3),'LineStyle','--')

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
hter=ternaryc(pMat(:,1),pMat(:,2),pMat(:,3));
set(hter,'marker','o','markerfacecolor','k','markersize',8,'markeredgecolor','k')


%% fig 2d: exp 1 model comparison
%% fig 3b: exp 2 main behavioral results
%% fig 4a: exp 2 modeling results

clear all; close all
fixedrisk = [];
colorMat = [1 0 0; 0 0 1; 0 0 0; 0 1 0];
nModelsPossible = 4; % how many total models there are
modelorderVec = [3 2 4 1]; % the order the models should be plotted
nModels = length(modelorderVec); % how many models you want to plot now
modelnameVec = {'max points','flexible','proportional','min error'};

nSubj = 11;

nPriorities = 3;
xlims = linspace(0,10,16); % for histograms
xdiff = diff(xlims(1:2));
filepath = 'fits/zuzanna/exp2/';

figure;
% top down, left right, bottom and top margins, left and right margins
ha = tight_subplot(3,nModels,{[.1 .1],[.03 .03 .03]},[.08 .03],[.05 .03]);
set(gcf,'Position',[28 504 800 500])


% load and plot predicted data for each model
for imodel = 1:nModels
    model = modelorderVec(imodel);
    
    % load data
    load([filepath 'modelpred_exp2_model' num2str(model) fixedrisk '.mat'],'preddata')
    
    % get model predictions
    for isubj = 1:nSubj
        for ipriority = 1:nPriorities
            % histogram of euclidean error
            simdatacounts = hist(preddata{isubj}{ipriority}(:,1),xlims);
            simerror{ipriority}(isubj,:) = simdatacounts./sum(simdatacounts);
            
            % histogram of circle size
            simdatacounts = hist(preddata{isubj}{ipriority}(:,2),xlims);
            simdiscsize{ipriority}(isubj,:)  = simdatacounts./sum(simdatacounts);
        end
    end
    
    meansimerror = cellfun(@mean,simerror,'UniformOutput',false);
    semsimerror = cellfun(@(x) std(x)./sqrt(size(x,1)),simerror,'UniformOutput',false);
    
    
    meansimdiscsize = cellfun(@mean,simdiscsize,'UniformOutput',false);
    semsimdiscsize = cellfun(@(x) std(x)./sqrt(size(x,1)),simdiscsize,'UniformOutput',false);
    
    
    % plot error (and circle size) model predictions
    for ipriority = 1:nPriorities
        
        axes(ha(imodel))
        fill( [0 xlims+xdiff fliplr(xlims)+xdiff 0],[0 meansimerror{ipriority}-semsimerror{ipriority}...
            fliplr(meansimerror{ipriority}+semsimerror{ipriority}) 0],colorMat(ipriority,:),'EdgeColor','none','FaceAlpha',0.4);
        hold on;
        title(modelnameVec{model})
        
        % circle size
        axes(ha(nModels+imodel))
        fill([0 xlims fliplr(xlims) 0]+xdiff,[0 meansimdiscsize{ipriority}-semsimdiscsize{ipriority}...
            fliplr(meansimdiscsize{ipriority}+semsimdiscsize{ipriority}) 0],colorMat(ipriority,:),'EdgeColor','none','FaceAlpha',0.4);
        hold on;
    end
    
    % ========== quantile correlation plot per subject ===========
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
    
    % ================ group correlation plot ====================
    
    meanmeanquantsimerror = cellfun(@nanmean,meanquantsimerror,'UniformOutput',false);
    meanmeanquantsimdiscsize = cellfun(@mean,meanquantsimdiscsize,'UniformOutput',false);
    semmeanquantsimdiscsize= cellfun(@(x) std(x)./sqrt(size(x,1)),meanquantsimdiscsize,'UniformOutput',false);
    
    axes(ha(2*nModels+imodel)); hold on
    for ipriority = 1:nPriorities
        plot_summaryfit(meanmeanquantsimerror{ipriority},[],[],meanmeanquantsimdiscsize{ipriority},...
            semmeanquantsimdiscsize{ipriority},[],colorMat(ipriority,:))
    end
    
end

% ============= PLOT REAL DATA ======================

% load data
load('exp2_cleandata.mat','data')
% load([filepath 'fits_model' num2str(model) '.mat'])


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

% =========== group plot =====================

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
    
    % ========== quantile correlation plot per subject ===========
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

%% fig 4b: exp 2 ternary plot

clear all
colorMat = [1 0 0; 0 0 1; 0 0 0];

% load data
load('fits/zuzanna/exp2/fits_model2.mat')
nSubj = size(ML_parameters,1);
pMat = ML_parameters(:,end-1:end);
pMat(:,3) = 1-sum(pMat,2);

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

% make patch showing monotonic area
x=0.5-c1*cos(pi/3)+c2/2;
y=0.866-c1*sin(pi/3)-c2*cot(pi/6)/2;
patch('Faces',[1 2 3],'Vertices',[x' y'],'FaceColor',grey,'FaceAlpha',0.3,'EdgeColor','none');

% change the color of the grid lines
set(hg(:,1),'color',axis1)
set(hg(:,2),'color',axis2)
set(hg(:,3),'color',axis3)

% make 0.6 0.3 0.1 lines noticeable
set(hg(3,1),'LineStyle','--')
set(hg(6,2),'LineStyle','--')
set(hg(1,3),'LineStyle','--')

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
hter=ternaryc(pMat(:,1),pMat(:,2),pMat(:,3));
set(hter,'marker','o','markerfacecolor','k','markersize',8,'markeredgecolor','k')


%% fig 4c: exp 2 model comparison


%% ====================================
%        SUPPLEMENTARY FIGURES
%  ====================================
%
% fig 1: how optimal allocation changes with Jbar for each model
% fig 2: how optimal allocation changes with experimental probe probability

%% fig 1:  how optimal allocation changes with Jbar for each model

clear all
% close all
exppriorityVec = [0.6 0.3 0.1];

% SET UP TERNARY PLOT
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
nSamps = 10;

Jbar = 2;
tau = 0.4;
alpha = 1;
beta = 1;
gamma = 1;

% Vec1 = linspace(0.01,3,nSamps); % tau
Vec2 = linspace(3,25,nSamps); % multiplier to use for Jbar


[pVec_MP, pVec_ME] = deal(nan(nSamps,3));
for isamp2 = 1:nSamps
    Jbar = Vec2(isamp2)*tau + 0.01
    
    MPtheta = [Jbar tau alpha beta];
    MEtheta = [MPtheta gamma];
    
    % maximizing points
    pVec_MP(isamp2,:) = calc_pVec_maxpoints(MPtheta,exppriorityVec);
    
    % minimizing error
    pVec_ME(isamp2,:) = calc_pVec_minerror(MEtheta,exppriorityVec);
end


sz = 40;
c = linspace(1,nSamps,nSamps);
% PLOT MAX POINTS MODEL
x=0.5-pVec_MP(:,1).*cos(pi/3)+pVec_MP(:,2)./2;
y=0.866-pVec_MP(:,1).*sin(pi/3)-pVec_MP(:,2).*cot(pi/6)/2;
plot(x,y,'k-')
colormap('autumn')
s = scatter(x,y,sz,c,'filled')

% PLOT MIN ERROR MODEL
x=0.5-pVec_ME(:,1).*cos(pi/3)+pVec_ME(:,2)./2;
y=0.866-pVec_ME(:,1).*sin(pi/3)-pVec_ME(:,2).*cot(pi/6)/2;
plot(x,y,'k-')
colormap('winter')
s2 = scatter(x,y,sz,c,'filled')


%% fig 2: how optimal allocation changes with experimental probe probability
% for minimizing error model

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
x = 0:0.2:1;
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
    exppriorityVec = [X(isamp) Y(isamp) Z(isamp)]
    
    pVec_ME(isamp,:) = calc_pVec_minerror(theta,exppriorityVec);
end

u = 0.5-pVec_ME(:,1).*cos(pi/3)+pVec_ME(:,2)./2;
v = 0.866-pVec_ME(:,1).*sin(pi/3)-pVec_ME(:,2).*cot(pi/6)/2;

u = u-x;
v = v-y;
quiver(x,y,u,v,Scale(0),'k')

%% how proportion allocation changes with probe probability
% for MP or ME model for each subject in experimet 2
% (not actually in supplementary...yet)

clear all
close all
clc

model = 'max_points';

% load MLEs for participants
switch model
    case 'max_points'
        load('fits/priority/exp2/fits_model1.mat')
    case 'min_error'
        load('fits/priotiy/exp2/fits_model2.mat')
end
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

% get x and y of actual values
x=0.5-X.*cos(pi/3)+Y./2;
y=0.866-X.*sin(pi/3)-Y.*cot(pi/6)/2;
nSamps = length(X);

figure;
for isubj = 1:nSubj;
    isubj
    
    % theta
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
    pVec_MP = nan(nSamps,3);
    for isamp = 1:nSamps;
        exppriorityVec = [X(isamp) Y(isamp) Z(isamp)];
        switch model
            case 'max_points'
                pVec_MP(isamp,:) = calc_pVec_maxpoints(theta,exppriorityVec);
            case 'min_error'
                pVec_MP(isamp,:) = calc_pVec_minerror(theta,exppriorityVec);
        end
    end
    
    u = 0.5-pVec_MP(:,1).*cos(pi/3)+pVec_MP(:,2)./2;
    v = 0.866-pVec_MP(:,1).*sin(pi/3)-pVec_MP(:,2).*cot(pi/6)/2;
    
    u = u-x;
    v = v-y;
    quiver(x,y,u,v,Scale(0),'k')
    
end