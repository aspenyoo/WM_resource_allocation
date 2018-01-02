


%% =================== EXP 1 ==========================
%           MAIN BEHAVIORAL RESULTS
% =====================================================

clear all
expnumber = 1;
load(['exp' num2str(expnumber) '_cleandata.mat'],'data')
filepath = ['fits/exp' num2str(expnumber) '/'];
nSubj = length(data);

priorityVec = [0.1 0.3 0.6];
nPriorities = 3;
colorMat = [0.5 0.5 0.5; 0.5 0.5 1; 1 0.5 0.5];
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

% get RT data
load(['exp' num2str(expnumber) '_zuzprocesseddata.mat'])
for isubj = 1:nSubj;
    for ipriority = 1:nPriorities;
        priority = priorityVec(ipriority);
        idx = (group_data(:,1) == isubj) & (group_data(:,2) == priority);
        SRTVec(isubj,ipriority) = mean(group_data(idx,3));
    end
end
% M and SEM
dataMSRT = nanmean(SRTVec);
dataSEMSRT = nanstd(SRTVec)./sqrt(nSubj);

dx = 0.2;
%  ========== plot ==========
subplot(1,2,1); hold on;
defaultplot;
axis([0.5 3.5 1 2.2])
set(gca,'YTick',1:.4:2.2,'YTickLabel',{1, '','', 2.2});
for ipriority = 1:nPriorities
    errorbar(ipriority,dataMerror(nPriorities+1-ipriority),dataSEMerror(nPriorities+1-ipriority),...
        'k','Marker',markerMat{ipriority},'LineWidth',1);
    errorb(ipriority,dataMerror(nPriorities+1-ipriority),dataSEMerror(nPriorities+1-ipriority),...
        'LineWidth',1);
end
xlabel('priority'); ylabel('error');
set(gca,'XTick',1:3,'XTickLabel',{'low','medium','high'})
title('error')

subplot(1,2,2); hold on;
defaultplot;
axis([0.5 3.5 430 530])
set(gca,'YTick',[430 480 530],'YTickLabel',{430 480 530});
xlabel('priority'); ylabel('error');
set(gca,'XTick',1:3,'XTickLabel',{'low','medium','high'})
title('SRT')
for ipriority = 1:nPriorities
    errorbar(ipriority,dataMSRT(ipriority),dataSEMSRT(ipriority),...
        'k','Marker',markerMat{ipriority},'LineWidth',1);
    errorb(ipriority,dataMSRT(ipriority),dataSEMSRT(ipriority),...
        'LineWidth',1);
end

%% =================== EXP 2 ==========================
%           MAIN BEHAVIORAL RESULTS
% =====================================================

clear all
expnumber = 2;
load(['exp' num2str(expnumber) '_cleandata.mat'],'data')
filepath = ['fits/exp' num2str(expnumber) '/'];
nSubj = length(data);

priorityVec = [0.1 0.3 0.6];
nPriorities = 3;
colorMat = [0 0 0; 0 0 1; 1 0 0];
markerMat = {'.','x','s'};

modelnameVec = {'optimal','free','fixed'}';

figure; hold on

% get error data
[dataerrorVec, datadiscsizeVec] = deal(nan(nSubj,nPriorities));
for ipriority = 1:nPriorities
    for isubj = 1:nSubj
        dataerrorVec(isubj,ipriority) = mean(data{isubj}{ipriority}(:,1));
        datadiscsizeVec(isubj,ipriority) = mean(data{isubj}{ipriority}(:,2));
    end
end
% M and SEM
dataMerror = mean(dataerrorVec);
dataSEMerror = std(dataerrorVec)./sqrt(nSubj);

dataMdiscsize = mean(datadiscsizeVec);
dataSEMdiscsize = std(datadiscsizeVec)./sqrt(nSubj);


dx = 0.2;
%  ========== plot ==========
subplot(1,2,1); hold on;
defaultplot;
xlabel('priority'); ylabel('error');
set(gca,'XTick',1:3,'XTickLabel',{'low','medium','high'})
title('error')
axis([0.5 3.5 0 2.4])
set(gca,'YTick',0:.4:2.4,'YTickLabel',0:.4:2.4);

for ipriority = 1:nPriorities
%     errorbar(ipriority,dataMerror(nPriorities+1-ipriority),dataSEMerror(nPriorities+1-ipriority),...
%         'LineColor',colorMat(ipriority,:),'LineWidth',1);
    errorb(ipriority,dataMerror(nPriorities+1-ipriority),dataSEMerror(nPriorities+1-ipriority),...
        'LineWidth',1,'Color',colorMat(ipriority,:));
end


subplot(1,2,2); hold on;
defaultplot;
axis([0.5 3.5 0 3.6])
set(gca,'YTick',0:0.6:3.6,'YTickLabel',0:0.6:3.6);
xlabel('priority'); ylabel('error');
set(gca,'XTick',1:3,'XTickLabel',{'low','medium','high'})
title('disc size')
for ipriority = 1:nPriorities
%     errorbar(ipriority,dataMdiscsize(nPriorities+1-ipriority),dataSEMdiscsize(nPriorities+1-ipriority),...
%         'k','Marker',markerMat{ipriority},'LineWidth',1);
    errorb(ipriority,dataMdiscsize(nPriorities+1-ipriority),dataSEMdiscsize(nPriorities+1-ipriority),...
        'LineWidth',1,'Color',colorMat(ipriority,:));
end

%% correlation plot exp 2

clear all
expnumber = 2;
imodel = 3;
fixedrisk = [];%'_fixedrisk';
loadpreddata = 1;
indvlplot = 0;

nPriorities = 3;
nTrials = 1e3*ones(1,3); % how many trials to simulate per priority
load(['exp' num2str(expnumber) '_cleandata.mat'],'data')
if (expnumber == 1)
    nSubj = 14;
else
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

% % histograms per subjects
% xlims = linspace(0,10,16);
%
% for isubj = 1:nSubj
%     if (indvlplot); figure; end;
%     for ipriority = 1:nPriorities
%
%         % histogram of euclidean error
%         datacounts = hist(data{isubj}{ipriority}(:,1),xlims);
%         simdatacounts = hist(preddata{isubj}{ipriority}(:,1),xlims);
%         error{ipriority}(isubj,:) = datacounts./sum(datacounts);
%         simerror{ipriority}(isubj,:) = simdatacounts./sum(simdatacounts);
%
%         if (indvlplot)
%             subplot(3,2,2*ipriority-1)
%             plot(xlims,error{ipriority}(isubj,:),'k')
%             hold on;
%             plot(xlims,simerror{ipriority}(isubj,:),'Color',aspencolors('booger'));
%             defaultplot
%             if ipriority == 1; title('euclidean error'); end
%         end
%
%         if (expnumber == 2)
%             % histogram of disc size
%             datacounts = hist(data{isubj}{ipriority}(:,2),xlims);
%             simdatacounts = hist(preddata{isubj}{ipriority}(:,2),xlims);
%             discsize{ipriority}(isubj,:) = datacounts./sum(datacounts);
%             simdiscsize{ipriority}(isubj,:)  = simdatacounts./sum(simdatacounts);
%
%             if (indvlplot)
%                 subplot(3,2,2*ipriority)
%                 plot(xlims,discsize{ipriority}(isubj,:),'k')
%                 hold on;
%                 plot(xlims,simdiscsize{ipriority}(isubj,:),'Color',aspencolors('booger'));
%                 defaultplot
%                 if ipriority == 1; title('disc size'); end
%             end
%         end
%     end
%     %     pause;
% end
%
% % =========== group plot =====================
%
% figure;
% colorMat = {'r','b','k'};
% if (expnumber == 2)
%     ha = tight_subplot(3,3,{[.03 .03],[.03 .07]},[.1 .01],[.1 .01]);
%     set(gcf,'Position',[28 504 800 500])
% else
%     ha = tight_subplot(1,3,.03,[.26 .05],[.11 .05]);
%     set(gcf,'Position',[28 504 800 236])
% end
%
%
% meanerror = cellfun(@mean,error,'UniformOutput',false);
% semerror = cellfun(@(x) std(x)./sqrt(size(x,1)),error,'UniformOutput',false);
% meansimerror = cellfun(@mean,simerror,'UniformOutput',false);
% semsimerror = cellfun(@(x) std(x)./sqrt(size(x,1)),simerror,'UniformOutput',false);
%
% if (expnumber == 2)
%     meandiscsize = cellfun(@mean,discsize,'UniformOutput',false);
%     semdiscsize = cellfun(@(x) std(x)./sqrt(size(x,1)),discsize,'UniformOutput',false);
%     meansimdiscsize = cellfun(@mean,simdiscsize,'UniformOutput',false);
%     semsimdiscsize = cellfun(@(x) std(x)./sqrt(size(x,1)),simdiscsize,'UniformOutput',false);
% end
%
% for ipriority = 1:nPriorities
%
%     if(expnumber == 2)
%         axes(ha(3*ipriority-2))
%     else
%         axes(ha(ipriority))
%     end
%     fill([xlims fliplr(xlims)],[meansimerror{ipriority}-semsimerror{ipriority}...
%         fliplr(meansimerror{ipriority}+semsimerror{ipriority})],colorMat{ipriority},'EdgeColor','none','FaceAlpha',0.4);
%     hold on;
%     errorbar(xlims,meanerror{ipriority},semerror{ipriority},'Color','k','LineStyle','none','LineWidth',1);
%     defaultplot
%     axis([0 10 0 0.4])
%
%     if expnumber == 2
%         %          axis([0 10 0 0.6])
%         if ipriority == 3
%             xlabel('error');
%         else
%             set(ha(3*ipriority-2),'XTickLabel','');
%         end
%         ylabel('proportion','FontSize',14);
%
%     else
%
%         xlabel('error','FontSize',16)
%         set(ha(ipriority),'YTick',[0 0.2 0.4],'FontSize',12);
%         if ipriority ~= 1
%             set(ha(ipriority),'YTickLabel','');
%         else
%             ylabel('proportion','FontSize',16)
%         end
%
%     end
%
%
%     if (expnumber == 2)
%         % discsize
%         axes(ha(3*ipriority-1))
%         fill([xlims fliplr(xlims)],[meansimdiscsize{ipriority}-semsimdiscsize{ipriority}...
%             fliplr(meansimdiscsize{ipriority}+semsimdiscsize{ipriority})],colorMat{ipriority},'EdgeColor','none','FaceAlpha',0.4);
%         hold on;
%         errorbar(xlims,meandiscsize{ipriority},semdiscsize{ipriority},'Color','k','LineStyle','none','LineWidth',1);
%         defaultplot
%         axis([0 10 0 0.6])
%         if ipriority == 3
%             xlabel('disc size','FontSize',14);
%         else
%             set(ha(3*ipriority-1),'XTickLabel','');
%         end
%         set(ha(3*ipriority-1),'YTickLabel','');
%
%     end
%
% end

% if (expnumber == 2)
% ========== quantile correlation plot per subject ===========
nQuants = 6;
for isubj = 1:nSubj
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
        end
        
    end
    
end

% ================ group plot ====================

meanmeanquanterror = cellfun(@mean,meanquanterror,'UniformOutput',false);
meanmeanquantdiscsize = cellfun(@mean,meanquantdiscsize,'UniformOutput',false);
semmeanquantdiscsize= cellfun(@(x) std(x)./sqrt(size(x,1)),meanquantdiscsize,'UniformOutput',false);

figure; hold on
colorMat = {'r','b','k'};
for ipriority = 1:nPriorities
    
    plot_summaryfit(meanmeanquanterror{ipriority},meanmeanquantdiscsize{ipriority},semmeanquantdiscsize{ipriority},...
        [],[],colorMat{ipriority})
    
    ylabel('disc size');
    axis([0 6 2 4])
    set(gca,'XTick',[0 3 6],'YTick',[2 3 4]);
    set(gca,'YTickLabel',[2 3 4])
    if (expnumber == 1)
        set(gca,'XTickLabel',[0 3 6])
        xlabel('error');
    end
end
set(gca,'XTickLabel',[0 3 6]);
xlabel('error');
% end

%% ==========   MAIN MODELING RESULTS  ==================

clear all
expnumber = 1;
imodel = 3;
fixedrisk = [];%'_fixedrisk';
loadpreddata = 1;
indvlplot = 0;

nPriorities = 3;
nTrials = 1e3*ones(1,3); % how many trials to simulate per priority
load(['exp' num2str(expnumber) '_cleandata.mat'],'data')
if (expnumber == 1)
    nSubj = 14;
else
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
    set(gcf,'Position',[28 504 800 500])
else
    ha = tight_subplot(1,3,.03,[.26 .05],[.11 .05]);
    set(gcf,'Position',[28 504 800 236])
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

if (expnumber == 2)
    % ========== quantile correlation plot per subject ===========
    nQuants = 6;
    for isubj = 1:nSubj
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
    set(ha(3*ipriority),'XTickLabel',[0 3 6]);
    xlabel('error');
end


%% =============== TERNARY PLOT ========================
clear all
% colorMat = [0.3 0.3 0.3; 0.3 0.3 1; 1 0.3 0.3];
colorMat = [1 0 0; 0 0 1; 0 0 0];

expnumber = 2;
imodel = 2;

% load data
filepath = ['fits/exp' num2str(expnumber) '/'];
load([filepath 'fits_model' num2str(imodel) '.mat'])
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
colorMap = aspencolors(nSubj,'qualitative');
for isubj = 1:nSubj
set(hter(isubj),'marker','o','markerfacecolor',colorMap(isubj,:),'markersize',8,'markeredgecolor',colorMap(isubj,:))
end
%% ====================================================
%          EXP 1: MAIN EFFECT AND MODELING RESULTS
% =====================================================

% clear all
% expnumber = 1;
% load(['exp' num2str(expnumber) '_cleandata.mat'],'data')
% filepath = ['fits/exp' num2str(expnumber) '/'];
% nSubj = length(data);
% 
% priorityVec = [0.1 0.3 0.6];
% nPriorities = 3;
% colorMat = [0.5 0.5 0.5; 0.5 0.5 1; 1 0.5 0.5];
% markerMat = {'.','x','s'};
% 
% modelnameVec = {'optimal','free','fixed'}';
% 
% figure; hold on
% 









clear all
expnumber = 1;
fixedrisk = [];%'_fixedrisk';
% colorMat = [aspencolors('dustyrose'); aspencolors('booger'); aspencolors('seacolored')];
colorMat = [1 0 0; 0 0 1; 0 0 0];
modelVec = [3 4 2];
nModels = length(modelVec);

if (expnumber == 1)
    nSubj = 14;
else
    nSubj = 11;
    
end
nPriorities = 3;
xlims = linspace(0,10,16); % for histograms
filepath = ['fits/exp' num2str(expnumber) fixedrisk '/'];

figure;
ha = tight_subplot(1,nModels+1,.03,[.26 .05],[.11 .05]);
set(gcf,'Position',[28 504 800 236])



% ----- REAL DATA STUFF -----

% load data
load(['exp' num2str(expnumber) '_cleandata.mat'],'data')

% histograms per subjects
for isubj = 1:nSubj
    for ipriority = 1:nPriorities
        
        % histogram of euclidean error
        datacounts = hist(data{isubj}{ipriority}(:,1),xlims);
        error{ipriority}(isubj,:) = datacounts./sum(datacounts);
        
        if (expnumber == 2)
            % histogram of disc size
            datacounts = hist(data{isubj}{ipriority}(:,2),xlims);
            discsize{ipriority}(isubj,:) = datacounts./sum(datacounts);
        end
    end
    
end

meanerror = cellfun(@mean,error,'UniformOutput',false);
semerror = cellfun(@(x) std(x)./sqrt(size(x,1)),error,'UniformOutput',false);


% ----- PLOT MAIN EFFECT ------

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

%  ========== plot ==========
axes(ha(1)); hold on;
defaultplot;
axis([0.5 3.5 1 2.2])
set(gca,'YTick',1:.4:2.2,'YTickLabel',{1, '','', 2.2});
for ipriority = 1:nPriorities
    errorb(ipriority,dataMerror(nPriorities+1-ipriority),dataSEMerror(nPriorities+1-ipriority),...
        'Color',colorMat(ipriority,:),'LineWidth',1);
end
xlabel('priority','FontSize',16); ylabel('error','FontSize',16);
set(ha(1),'XTick',1:3,'XTickLabel',{'low','medium','high'})


% ------- PRED DATA STUFF --------

% load and plot predicted data for each model
for imodel = 1:nModels
    modelnum = modelVec(imodel);
    
    % load data
    load([filepath 'modelpred_exp' num2str(expnumber) '_model' num2str(modelnum) fixedrisk '.mat'],'preddata')
    
    for isubj = 1:nSubj
        for ipriority = 1:nPriorities
            
            % histogram of euclidean error
            simdatacounts = hist(preddata{isubj}{ipriority}(:,1),xlims);
            simerror{ipriority}(isubj,:) = simdatacounts./sum(simdatacounts);
            
            if (expnumber == 2)
                % histogram of disc size
                simdatacounts = hist(preddata{isubj}{ipriority}(:,2),xlims);
                simdiscsize{ipriority}(isubj,:)  = simdatacounts./sum(simdatacounts);
                
            end
        end
    end
    
    meansimerror = cellfun(@mean,simerror,'UniformOutput',false);
    semsimerror = cellfun(@(x) std(x)./sqrt(size(x,1)),simerror,'UniformOutput',false);
    
    if (expnumber == 2)
        meansimdiscsize = cellfun(@mean,simdiscsize,'UniformOutput',false);
        semsimdiscsize = cellfun(@(x) std(x)./sqrt(size(x,1)),simdiscsize,'UniformOutput',false);
    end
    
    axes(ha(imodel+1)); hold on;
    for ipriority = 1:nPriorities
        fill([xlims fliplr(xlims)],[meansimerror{ipriority}-semsimerror{ipriority}...
            fliplr(meansimerror{ipriority}+semsimerror{ipriority})],colorMat(ipriority,:),'EdgeColor','none','FaceAlpha',0.4);
        
        errorbar(xlims,meanerror{ipriority},semerror{ipriority},'Color',colorMat(ipriority,:),'LineStyle','none','LineWidth',1);
    end
    defaultplot
    axis([0 10 0 0.4])
    xlabel('error','FontSize',16)
    set(ha(imodel+1),'YTick',[0 0.2 0.4],'YTickLabel',[0 0.2 0.4],'XTick',[0 5 10],'XTickLabel',[0 5 10]);
    if imodel ~= 1
        set(ha(imodel+1),'YTickLabel','');
    else
        ylabel('proportion','FontSize',16)
    end
        
    
end

%% ====================================================
%          MODELING RESULTS: COMBINED PRIORITIES
% =====================================================

clear all; close all
expnumber = 1;
fixedrisk = [];  
colorMat = [1 0 0; 0 0 1; 0 0 0; 0 1 0];
nModelsPossible = 4; % how many total models there are
modelorderVec = [3 2 4]; % the order the models should be plotted
nModels = length(modelorderVec); % how many models you want to plot now
modelnameVec = {'max points','flexible','proportional','min error'};

if (expnumber == 1)
    nSubj = 14;
else
    nSubj = 11;
end
nPriorities = 3;
xlims = linspace(0,10,16); % for histograms
xdiff = diff(xlims(1:2));
filepath = ['fits/exp' num2str(expnumber) fixedrisk '/'];

figure;
if (expnumber == 2)
    % top down, left right, bottom and top margins, left and right margins
    ha = tight_subplot(3,nModels,{[.1 .1],[.03 .03 .03]},[.08 .03],[.05 .03]);
    set(gcf,'Position',[28 504 800 500])
else
    ha = tight_subplot(1,nModels,[.03 .03],[.2 .07],[.06 .01]);
    set(gcf,'Position',[28 504 800 236])
end

% load and plot predicted data for each model
for imodel = 1:nModels
    model = modelorderVec(imodel);
    
    % load data
    load([filepath 'modelpred_exp' num2str(expnumber) '_model' num2str(model) fixedrisk '.mat'],'preddata')
    
    % get model predictions
    for isubj = 1:nSubj
        for ipriority = 1:nPriorities
            % histogram of euclidean error
            simdatacounts = hist(preddata{isubj}{ipriority}(:,1),xlims);
            simerror{ipriority}(isubj,:) = simdatacounts./sum(simdatacounts);
            if (expnumber == 2)
                % histogram of disc size
                simdatacounts = hist(preddata{isubj}{ipriority}(:,2),xlims);
                simdiscsize{ipriority}(isubj,:)  = simdatacounts./sum(simdatacounts);
            end
        end
    end
    
    meansimerror = cellfun(@mean,simerror,'UniformOutput',false);
    semsimerror = cellfun(@(x) std(x)./sqrt(size(x,1)),simerror,'UniformOutput',false);
    
    if (expnumber == 2)
        meansimdiscsize = cellfun(@mean,simdiscsize,'UniformOutput',false);
        semsimdiscsize = cellfun(@(x) std(x)./sqrt(size(x,1)),simdiscsize,'UniformOutput',false);
    end
    
    % plot error (and disc size) model predictions
    for ipriority = 1:nPriorities
        
            axes(ha(imodel))
        fill( [0 xlims+xdiff fliplr(xlims)+xdiff 0],[0 meansimerror{ipriority}-semsimerror{ipriority}...
            fliplr(meansimerror{ipriority}+semsimerror{ipriority}) 0],colorMat(ipriority,:),'EdgeColor','none','FaceAlpha',0.4);
        hold on;
        title(modelnameVec{model})
        
        if (expnumber == 2)
            % discsize
            axes(ha(nModels+imodel))
            fill([0 xlims fliplr(xlims) 0]+xdiff,[0 meansimdiscsize{ipriority}-semsimdiscsize{ipriority}...
                fliplr(meansimdiscsize{ipriority}+semsimdiscsize{ipriority}) 0],colorMat(ipriority,:),'EdgeColor','none','FaceAlpha',0.4);
            hold on;
        end
        
    end
    
    if (expnumber == 2)
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
    
end

% ============= PLOT REAL DATA ======================

% load data
load(['exp' num2str(expnumber) '_cleandata.mat'],'data')
% load([filepath 'fits_model' num2str(model) '.mat'])


% histograms per subjects
for isubj = 1:nSubj
    for ipriority = 1:nPriorities
        
        % histogram of euclidean error
        datacounts = hist(data{isubj}{ipriority}(:,1),xlims);
        error{ipriority}(isubj,:) = datacounts./sum(datacounts);
        
        if (expnumber == 2)
            % histogram of disc size
            datacounts = hist(data{isubj}{ipriority}(:,2),xlims);
            discsize{ipriority}(isubj,:) = datacounts./sum(datacounts);
        end
    end
    
end

% =========== group plot =====================

for imodel = 1:nModels

meanerror = cellfun(@mean,error,'UniformOutput',false);
semerror = cellfun(@(x) std(x)./sqrt(size(x,1)),error,'UniformOutput',false);

if (expnumber == 2)
    meandiscsize = cellfun(@mean,discsize,'UniformOutput',false);
    semdiscsize = cellfun(@(x) std(x)./sqrt(size(x,1)),discsize,'UniformOutput',false);
end

for ipriority = 1:nPriorities
    
    
    % =========== ERROR ============
    
        axes(ha(imodel))
    hold on;
    errorbar(xlims+xdiff,meanerror{ipriority},semerror{ipriority},'Color',colorMat(ipriority,:),'LineStyle', 'none','LineWidth',1);
    defaultplot
    axis([0 10 0 0.4])
    
    if expnumber == 2
        xlabel('error');
        set(ha(imodel),'XTick',0:10,'YTick',0:0.1:0.4,...
            'YTickLabel','','FontSize',10);
        if (imodel == 1)
        ylabel('proportion','FontSize',14);
        set(ha(imodel),'YTickLabel',0:0.1:0.4)
        end
    else
        xlabel('error','FontSize',16)
        set(ha(imodel),'YTick',0:0.1:0.4,'XTick',0:10);
        
        if imodel == 2
            ylabel('proportion','FontSize',16)
        else
            set(ha(imodel),'YTickLabel','');
        end
    end
    
    % ========= DISC SIZE ===========
    if (expnumber == 2)
        axes(ha(nModels+imodel))
        hold on;
        errorbar(xlims+xdiff,meandiscsize{ipriority},semdiscsize{ipriority},'Color',colorMat(ipriority,:),'LineStyle','none','LineWidth',1);
        defaultplot
        axis([0 10 0 0.4])
        set(ha(nModels+imodel),'XTick',0:10,'YTick',0:0.1:0.4,...
            'FontSize',10,'YTickLabel','');
        xlabel('disc size','FontSize',10);
        if (imodel == 1); 
            ylabel('proportion','FontSize',14); 
            set(ha(nModels+imodel),'YTickLabel',0:0.1:0.4);
        end
    end
end

if (expnumber == 2)
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
%         plot_summaryfit(meanmeanquanterror{ipriority},...
%             [],[],colorMat(ipriority,:))
    end
    
    if (imodel == 1); 
        ylabel('disc size','FontSize',14); 
    end
    set(ha(2*nModels+imodel),'XTick',0:6,'YTick',0:6,'XTickLabel',0:6,'YTickLabel',0:6);
   
    if (expnumber == 1)
%         set(ha(imodel),'XTickLabel',[0 3 6])
%         xlabel('error');
    end
    
    xlabel('error');
end
end

%% ========================================================================
%      MODELING RESULTS: COMBINED PRIORITIES CUMULATIVE DISTRIBUTIONS
% =========================================================================

clear all; close all
expnumber = 2;
fixedrisk = [];  
colorMat = [1 0 0; 0 0 1; 0 0 0; 0 1 0];
nModelsPossible = 4; % how many total models there are
modelorderVec = [3 2 4 1]; % the order the models should be plotted
nModels = length(modelorderVec); % how many models you want to plot now
modelnameVec = {'max points','flexible','proportional','min error'};

if (expnumber == 1)
    nSubj = 14;
else
    nSubj = 11;
end
nPriorities = 3;
xlims = linspace(0,10,100); % for histograms
dataxlims = linspace(0,10,16);
filepath = ['fits/exp' num2str(expnumber) fixedrisk '/'];

figure;
if (expnumber == 2)
    % top down, left right, bottom and top margins, left and right margins
    ha = tight_subplot(3,nModels,{[.1 .1],[.03 .03 .03]},[.08 .03],[.05 .03]);
    set(gcf,'Position',[28 504 800 500])
else
    ha = tight_subplot(1,nModels,[.03 .03],[.2 .07],[.06 .01]);
    set(gcf,'Position',[28 504 800 236])
end

% load and plot predicted data for each model
for imodel = 1:nModels
    model = modelorderVec(imodel);
    
    % load data
    load([filepath 'modelpred_exp' num2str(expnumber) '_model' num2str(model) fixedrisk '.mat'],'preddata')
    
    % get model predictions
    for isubj = 1:nSubj
        for ipriority = 1:nPriorities
            % histogram of euclidean error
            simdatacounts = hist(preddata{isubj}{ipriority}(:,1),xlims);
            simerror{ipriority}(isubj,:) = cumsum(simdatacounts)./sum(simdatacounts);
            if (expnumber == 2)
                % histogram of disc size
                simdatacounts = hist(preddata{isubj}{ipriority}(:,2),xlims);
                simdiscsize{ipriority}(isubj,:)  = cumsum(simdatacounts)./sum(simdatacounts);
            end
        end
    end
    
    meansimerror = cellfun(@mean,simerror,'UniformOutput',false);
    semsimerror = cellfun(@(x) std(x)./sqrt(size(x,1)),simerror,'UniformOutput',false);
    
    if (expnumber == 2)
        meansimdiscsize = cellfun(@mean,simdiscsize,'UniformOutput',false);
        semsimdiscsize = cellfun(@(x) std(x)./sqrt(size(x,1)),simdiscsize,'UniformOutput',false);
    end
    
    % plot error (and disc size) model predictions
    for ipriority = 1:nPriorities
        
            axes(ha(imodel))
        fill([xlims fliplr(xlims)],[meansimerror{ipriority}-semsimerror{ipriority}...
            fliplr(meansimerror{ipriority}+semsimerror{ipriority})],colorMat(ipriority,:),'EdgeColor','none','FaceAlpha',0.4);
        hold on;
        title(modelnameVec{model})
        
        if (expnumber == 2)
            % discsize
            axes(ha(nModels+imodel))
            fill([xlims fliplr(xlims)],[meansimdiscsize{ipriority}-semsimdiscsize{ipriority}...
                fliplr(meansimdiscsize{ipriority}+semsimdiscsize{ipriority})],colorMat(ipriority,:),'EdgeColor','none','FaceAlpha',0.4);
            hold on;
        end
        
    end
    
    if (expnumber == 2)
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
    
end

% ============= PLOT REAL DATA ======================

% load data
load(['exp' num2str(expnumber) '_cleandata.mat'],'data')
% load([filepath 'fits_model' num2str(model) '.mat'])


% histograms per subjects
for isubj = 1:nSubj
    for ipriority = 1:nPriorities
        
        % histogram of euclidean error
        datacounts = hist(data{isubj}{ipriority}(:,1),dataxlims);
        error{ipriority}(isubj,:) = cumsum(datacounts)./sum(datacounts);
        
        if (expnumber == 2)
            % histogram of disc size
            datacounts = hist(data{isubj}{ipriority}(:,2),dataxlims);
            discsize{ipriority}(isubj,:) = cumsum(datacounts)./sum(datacounts);
        end
    end
    
end

% =========== group plot =====================

for imodel = 1:nModels

meanerror = cellfun(@mean,error,'UniformOutput',false);
semerror = cellfun(@(x) std(x)./sqrt(size(x,1)),error,'UniformOutput',false);

if (expnumber == 2)
    meandiscsize = cellfun(@mean,discsize,'UniformOutput',false);
    semdiscsize = cellfun(@(x) std(x)./sqrt(size(x,1)),discsize,'UniformOutput',false);
end

for ipriority = 1:nPriorities
    
    
    % =========== ERROR ============
    
        axes(ha(imodel))
    hold on;
%     plot(xlims,meanerror{ipriority},'--','Color',colorMat(ipriority,:));
%     plot(xlims,meanerror{ipriority}-semerror{ipriority},'Color',colorMat(ipriority,:));
%     plot(xlims,meanerror{ipriority}+semerror{ipriority},'Color',colorMat(ipriority,:));
    errorbar(dataxlims,meanerror{ipriority},semerror{ipriority},'Color',colorMat(ipriority,:),'LineWidth',1);
    defaultplot
    axis([0 10 0 1])
    
    if expnumber == 2
        xlabel('error');
        set(ha(imodel),'XTick',[0 5 10],'YTick',[0 0.5 1],...
            'YTickLabel','','FontSize',10);
        if (imodel == 1)
        ylabel('proportion','FontSize',14);
        set(ha(imodel),'YTickLabel',[0 0.5 1])
        end
    else
        xlabel('error','FontSize',16)
        set(ha(imodel),'YTick',[0 0.5 1],'XTick',[0 5 10]);
        
        if imodel == 2
            ylabel('proportion','FontSize',16)
        else
            set(ha(imodel),'YTickLabel','');
        end
    end
    
    % ========= DISC SIZE ===========
    if (expnumber == 2)
        axes(ha(nModels+imodel))
        hold on;
%         plot(xlims,meandiscsize{ipriority},'--','Color',colorMat(ipriority,:));
%         plot(xlims,meandiscsize{ipriority}-semdiscsize{ipriority},'Color',colorMat(ipriority,:));
%         plot(xlims,meandiscsize{ipriority}+semdiscsize{ipriority},'Color',colorMat(ipriority,:));
        errorbar(dataxlims,meandiscsize{ipriority},semdiscsize{ipriority},'Color',colorMat(ipriority,:),'LineWidth',1);
        defaultplot
        axis([0 10 0 1])
        set(ha(nModels+imodel),'XTick',[0 5 10],'YTick',[0 0.5 1],...
            'FontSize',10,'YTickLabel','');
        xlabel('disc size','FontSize',10);
        if (imodel == 1); 
            ylabel('proportion','FontSize',14); 
            set(ha(nModels+imodel),'YTickLabel',[0 0.5 1]);
        end
    end
end

if (expnumber == 2)
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
    for ipriority = 1:nPriorities
        hold on
        plot_summaryfit(meanmeanquanterror{ipriority},meanmeanquantdiscsize{ipriority},semmeanquantdiscsize{ipriority},...
            [],[],colorMat(ipriority,:))
    end
    
    if (imodel == 1); 
        ylabel('disc size','FontSize',14); 
    end
    axis([0 6 2 6])
    set(ha(2*nModels+imodel),'XTick',[0 3 6],'YTick',[2 4 6],...
        'XTickLabel',[0 3 6],'YTickLabel',[2 4 6]);
   
    if (expnumber == 1)
%         set(ha(imodel),'XTickLabel',[0 3 6])
%         xlabel('error');
    end
    
    xlabel('error');
end
end

%% ================================================================
%          MODELING RESULTS: COMBINED PRIORIOTIES MODELS AS ROWS
% =================================================================

clear all; close all
expnumber = 2;
fixedrisk = [];  %'_fixedrisk';
% colorMat = [aspencolors('dustyrose'); aspencolors('booger'); aspencolors('seacolored')];
colorMat = [1 0 0; 0 0 1; 0 0 0; 0 1 0];
nModelsPossible = 4; % how many total models there are
modelorderVec = [3 2]; % the order the models should be plotted
nModels = length(modelorderVec); % how many models you want to plot now
modelnameVec = {'max points','free','proportional','min error'};

if (expnumber == 1)
    nSubj = 14;
else
    nSubj = 11;
end
nPriorities = 3;
xlims = linspace(0,10,16); % for histograms
filepath = ['fits/exp' num2str(expnumber) fixedrisk '/'];

figure;
if (expnumber == 2)
    % top down, left right, bottom and top margins, left and right margins
    ha = tight_subplot(3,nModels,{[.03 .03],.05*ones(1,nModels-1)},[.08 .03],[.05 .03]);
    set(gcf,'Position',[28 504 800 500])
else
    ha = tight_subplot(1,nModels,[.03 .03],[.2 .07],[.06 .01]);
    set(gcf,'Position',[28 504 800 236])
end

% load and plot predicted data for each model
for imodel = 1:nModels
    model = modelorderVec(imodel);
    
    % load data
    load([filepath 'modelpred_exp' num2str(expnumber) '_model' num2str(model) fixedrisk '.mat'],'preddata')
    
    % get model predictions
    for isubj = 1:nSubj
        for ipriority = 1:nPriorities
            % histogram of euclidean error
            simdatacounts = hist(preddata{isubj}{ipriority}(:,1),xlims);
            simerror{ipriority}(isubj,:) = simdatacounts./sum(simdatacounts);
            if (expnumber == 2)
                % histogram of disc size
                simdatacounts = hist(preddata{isubj}{ipriority}(:,2),xlims);
                simdiscsize{ipriority}(isubj,:)  = simdatacounts./sum(simdatacounts);
            end
        end
    end
    
    meansimerror = cellfun(@mean,simerror,'UniformOutput',false);
    semsimerror = cellfun(@(x) std(x)./sqrt(size(x,1)),simerror,'UniformOutput',false);
    
    if (expnumber == 2)
        meansimdiscsize = cellfun(@mean,simdiscsize,'UniformOutput',false);
        semsimdiscsize = cellfun(@(x) std(x)./sqrt(size(x,1)),simdiscsize,'UniformOutput',false);
    end
    
    % plot error (and disc size) model predictions
    for ipriority = 1:nPriorities
        
            axes(ha(nModels*(imodel-1)+1))
        fill([xlims fliplr(xlims)],[meansimerror{ipriority}-semsimerror{ipriority}...
            fliplr(meansimerror{ipriority}+semsimerror{ipriority})],colorMat(ipriority,:),'EdgeColor','none','FaceAlpha',0.4);
        hold on;
        title(modelnameVec{model})
        
        if (expnumber == 2)
            % discsize
            axes(ha(nModels*(imodel-1)+2))
            fill([xlims fliplr(xlims)],[meansimdiscsize{ipriority}-semsimdiscsize{ipriority}...
                fliplr(meansimdiscsize{ipriority}+semsimdiscsize{ipriority})],colorMat(ipriority,:),'EdgeColor','none','FaceAlpha',0.4);
            hold on;
        end
        
    end
    
    if (expnumber == 2)
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
        
        axes(ha(nModels*(imodel-1)+3)); hold on
        for ipriority = 1:nPriorities
            plot_summaryfit(meanmeanquantsimerror{ipriority},[],[],meanmeanquantsimdiscsize{ipriority},...
                semmeanquantsimdiscsize{ipriority},[],colorMat(ipriority,:))
        end
    end
    
end

% ============= PLOT REAL DATA ======================

% load data
load(['exp' num2str(expnumber) '_cleandata.mat'],'data')
% load([filepath 'fits_model' num2str(model) '.mat'])


% histograms per subjects
for isubj = 1:nSubj
    for ipriority = 1:nPriorities
        
        % histogram of euclidean error
        datacounts = hist(data{isubj}{ipriority}(:,1),xlims);
        error{ipriority}(isubj,:) = datacounts./sum(datacounts);
        
        if (expnumber == 2)
            % histogram of disc size
            datacounts = hist(data{isubj}{ipriority}(:,2),xlims);
            discsize{ipriority}(isubj,:) = datacounts./sum(datacounts);
        end
    end
    
end

% =========== group plot =====================

for imodel = 1:nModels

meanerror = cellfun(@mean,error,'UniformOutput',false);
semerror = cellfun(@(x) std(x)./sqrt(size(x,1)),error,'UniformOutput',false);

if (expnumber == 2)
    meandiscsize = cellfun(@mean,discsize,'UniformOutput',false);
    semdiscsize = cellfun(@(x) std(x)./sqrt(size(x,1)),discsize,'UniformOutput',false);
end

for ipriority = 1:nPriorities
    
    
    % =========== ERROR ============
    
    axes(ha(nModels*(imodel-1)+1))
    hold on;
    errorbar(xlims,meanerror{ipriority},semerror{ipriority},'Color',colorMat(ipriority,:),'LineStyle','none','LineWidth',1);
    defaultplot
    axis([0 10 0 0.4])
    
    if expnumber == 2
        ylabel('proportion','FontSize',16);
        set(ha(nModels*(imodel-1)+1),'YTick',[0 0.2 0.4],...
            'XTickLabel','','XTick',[0 5 10])
        
        if (imodel == nModels)
        xlabel('error');
        set(ha(nModels*(imodel-1)+1),'YTick',[0 0.2 0.4],'XTickLabel',[0 5 10]);
        end
    else
        xlabel('error','FontSize',16)
        set(ha(nModels*(imodel-1)+1),'YTick',[0 0.2 0.4],'FontSize',12);
        
        if imodel == 2
            ylabel('proportion','FontSize',16)
        else
            set(ha(nModels*(imodel-1)+1),'YTickLabel','');
        end
    end
    
    % ========= DISC SIZE ===========
    if (expnumber == 2)
        axes(ha(nModels*(imodel-1)+2))
        hold on;
        errorbar(xlims,meandiscsize{ipriority},semdiscsize{ipriority},'Color',colorMat(ipriority,:),'LineStyle','none','LineWidth',1);
        defaultplot
        axis([0 10 0 0.4])
        ylabel('proportion','FontSize',16);
        set(ha(nModels*(imodel-1)+2),'YTick',[0 0.2 0.4],...
            'XTickLabel','','XTick',[0 5 10])
        if (imodel == nModels)
        xlabel('error');
        set(ha(nModels*(imodel-1)+2),'YTick',[0 0.2 0.4],'XTickLabel',[0 5 10]);
        end
    end
end

if (expnumber == 2)
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
    
    axes(ha(nModels*(imodel-1)+3))
    for ipriority = 1:nPriorities
        hold on
        plot_summaryfit(meanmeanquanterror{ipriority},meanmeanquantdiscsize{ipriority},semmeanquantdiscsize{ipriority},...
            [],[],colorMat(ipriority,:))
    end
    
    
    axis([0 6 2 6])
    ylabel('disc size','FontSize',16); 
    set(ha(nModels*(imodel-1)+3),'XTick',[0 3 6],'YTick',[2 4 6],...
        'XTickLabel','','YTickLabel',[2 4 6]);
    if (imodel == nModels); 
        xlabel('error');
        set(ha(nModels*(imodel-1)+3),'XTickLabel',[0 3 6])
    end
    
    
   
    if (expnumber == 1)
%         set(ha(imodel),'XTickLabel',[0 3 6])
%         xlabel('error');
    end
    
    
end
end

%% JUST DATA plot


clear all; close all
expnumber = 1;
fixedrisk = [];  %'_fixedrisk';
% colorMat = [aspencolors('dustyrose'); aspencolors('booger'); aspencolors('seacolored')];
colorMat = [1 0 0; 0 0 1; 0 0 0; 0 1 0];

if (expnumber == 1)
    nSubj = 14;
else
    nSubj = 11;
end
nPriorities = 3;
xlims = linspace(0,10,16); % for histograms
filepath = ['fits/exp' num2str(expnumber) fixedrisk '/'];

if (expnumber == 2)
    figure;
    % top down, left right, bottom and top margins, left and right margins
    ha = tight_subplot(1,3,{[.1],[.03 .03]},[.08 .03],[.05 .03]);
    set(gcf,'Position',[28 504 800 250])
else
    ha = tight_subplot(1,1,{[.1],[.03]},[.1 .04],[.15 .03]);
    set(gcf,'Position',[28 504 300 236])
end

% ============= PLOT REAL DATA ======================

% load data
load(['exp' num2str(expnumber) '_cleandata.mat'],'data')
% load([filepath 'fits_model' num2str(model) '.mat'])


% histograms per subjects
for isubj = 1:nSubj
    for ipriority = 1:nPriorities
        
        % histogram of euclidean error
        datacounts = hist(data{isubj}{ipriority}(:,1),xlims);
        error{ipriority}(isubj,:) = datacounts./sum(datacounts);
        
        if (expnumber == 2)
            % histogram of disc size
            datacounts = hist(data{isubj}{ipriority}(:,2),xlims);
            discsize{ipriority}(isubj,:) = datacounts./sum(datacounts);
        end
    end
    
end

% =========== group plot =====================

meanerror = cellfun(@mean,error,'UniformOutput',false);
semerror = cellfun(@(x) std(x)./sqrt(size(x,1)),error,'UniformOutput',false);
meanmeanerror = cellfun(@(x) sum(xlims.*x),meanerror,'UniformOutput',false);

if (expnumber == 2)
    meandiscsize = cellfun(@mean,discsize,'UniformOutput',false);
    semdiscsize = cellfun(@(x) std(x)./sqrt(size(x,1)),discsize,'UniformOutput',false);
    meanmeandiscsize = cellfun(@(x) sum(xlims.*x),meandiscsize,'UniformOutput',false);
end

for ipriority = 1:nPriorities
    
    % =========== ERROR ============
    
    axes(ha(1))
    hold on;
    errorbar(xlims,meanerror{ipriority},semerror{ipriority},'Color',colorMat(ipriority,:),'LineWidth',1);
    plot([meanmeanerror{ipriority} meanmeanerror{ipriority}],[0 0.4],'Color',colorMat(ipriority,:));
    defaultplot
    axis([0 10 0 0.4])
    
    ylabel('proportion','FontSize',16);
        set(ha(1),'YTick',[0 0.2 0.4],'YTickLabel',[0 0.2 0.4],...
            'XTick',[0 5 10],'XTickLabel',[0 5 10])
        xlabel('error','FontSize',16);
    
    % ========= DISC SIZE ===========
    if (expnumber == 2)
        axes(ha(2))
        hold on;
        errorbar(xlims,meandiscsize{ipriority},semdiscsize{ipriority},'Color',colorMat(ipriority,:),'LineWidth',1);
        plot([meanmeandiscsize{ipriority} meanmeandiscsize{ipriority}],[0 0.4],'Color',colorMat(ipriority,:));
        defaultplot
        axis([0 10 0 0.4])
        ylabel('proportion','FontSize',16);
        set(ha(2),'YTick',[0 0.2 0.4],'YTickLabel',[0 0.2 0.4],...
            'XTick',[0 5 10],'XTickLabel',[0 5 10])
        xlabel('error','FontSize',16);
    end
end

if (expnumber == 2)
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
    
    axes(ha(3))
    for ipriority = 1:nPriorities
        hold on
        errorbar(meanmeanquanterror{ipriority},meanmeanquantdiscsize{ipriority},semmeanquantdiscsize{ipriority},'Color',colorMat(ipriority,:),'LineWidth',1);
    end
    
    defaultplot
    axis([0 6 2 6])
    ylabel('disc size','FontSize',16); 
    set(ha(3),'XTick',[0 3 6],'YTick',[2 4 6],...
        'YTickLabel',[2 4 6],'XTickLabel',[0 5 10]);
        xlabel('error','FontSize',16);
    
    
   
    if (expnumber == 1)
%         set(ha(imodel),'XTickLabel',[0 3 6])
%         xlabel('error');
    end
    
    
end

%% schematic of VP model

% stimulus information
angles = [50 150 240 300]; % locations in degrees
variances = [.5 1 2]; % variance in dva

mean_x = 10.*cosd(angles);
mean_y = 10.*sind(angles);

nGrids = 100;
x = linspace(-12,12,nGrids);
[xx, yy] = meshgrid(x,x);

map = zeros(nGrids^2,1);
for istim = 1:(length(angles)-1);
    MU = [mean_x(istim) mean_y(istim)];
    SIGMA = [variances(istim) 0; 0 variances(istim)];
    
    pdf = mvnpdf([xx(:) yy(:)], MU, SIGMA);
    map = map + pdf(:);
end

contour3(reshape(map,nGrids,nGrids),50)
% plot3(xx(:),yy(:),map)

%% figure of how resource allocation changes as a function of Jbar_total
clear all
close all

nJbars = 5;
tau = 0.1;
JbarVec = exp(linspace(log(3*tau+0.01),log(15),nJbars));
alpha = 1; 
beta = 2; 

[pVec_MP, pVec_ME] = deal(nan(nJbars,3));
for iJbar = 1:nJbars;
    iJbar
    
    theta = [JbarVec(iJbar) tau alpha beta];
    
    % maxmimizing points
    pVec_MP(iJbar,:) = calc_optimal_pVec(theta);
    
    % minimizing error
    pVec_ME(iJbar,:) = calc_pVec_minerror(theta);
    
end


% =============== TERNARY PLOT ========================
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
set(h,'marker','.','markerfacecolor',[1 0.5 0],'markersize',24,'markeredgecolor',[1 0.5 0])

% plot MP 
% colormap(aspencolors(9,'yellowgreen'));
hter = ternaryc(pVec_MP(:,1),pVec_MP(:,2),pVec_MP(:,3),1:nJbars,'*');
hold on
for i = 1:(nJbars-1);
    plot([hter(i).XData hter(i+1).XData], [hter(i).YData hter(i+1).YData],'k-')
end

% plot ME 
hter = ternaryc(pVec_ME(:,1),pVec_ME(:,2),pVec_ME(:,3),1:nJbars,'o');
hold on
for i = 1:(nJbars-1);
    plot([hter(i).XData hter(i+1).XData], [hter(i).YData hter(i+1).YData],'k-')
end

%% PERMUTATION TEST FOR STIMULUS: EXP 2
clear all
close all

nBins = 24;
expnumber = 2; 

switch expnumber
    case 1
        load('exp1_zuzprocesseddata.mat')
        dataa = sqrt((group_data(:,13) - group_data(:,11)).^2 +...
            (group_data(:,14) - group_data(:,12)).^2);
        anglee = atand(group_data(:,14)./group_data(:,13));

        idx = (group_data(:,13) < 0) & (group_data(:,14) > 0); % second quadrant
        anglee(idx) = anglee(idx) + 180;
        
        idx = (group_data(:,13) < 0) & (group_data(:,14) < 0); % third quadrant
        anglee(idx) = anglee(idx) + 180;
        
        idx = (group_data(:,13) > 0) & (group_data(:,14) < 0); % third quadrant
        anglee(idx) = anglee(idx) + 360;

        group_data = [group_data(:,1) anglee dataa];
    case 2
        load('exp2_zuzprocesseddata.mat')
        group_data = group_data(:,[1 6 7 13]);
        subjVec = 4:14;
        dataa = sqrt((10.*cosd(group_data(:,4)) - group_data(:,2)).^2 +...
            (10.*sind(group_data(:,4)) - group_data(:,3)).^2);
        group_data = [group_data(:,[1 4]) dataa];
end
nSubj = length(subjVec);

group_data(any(isnan(group_data), 2), :) = [];

[angleMat, errorMat, permerrorMat, permangleMat] = deal(nan(nSubj,nBins));

for isubj = 1:nSubj;
    subjnum = subjVec(isubj);
    
    idx = group_data(:,1) == subjnum;
    nTrials = sum(idx);
    trialss = round(linspace(1,nTrials,nBins+1));
    
    dataa = group_data(idx,2:3);
    dataa = sortrows(dataa,1);
    
    permuteidx = randperm(nTrials);
    permutedataa = [dataa(:,1) dataa(permuteidx,2)];
    
    for ibin = 1:nBins;
        bindata = dataa(trialss(ibin):trialss(ibin+1),:);
        permdata = permutedataa(trialss(ibin):trialss(ibin+1),:);
        
        angleMat(isubj,ibin) = mean(bindata(:,1));
        errorMat(isubj,ibin) = mean(bindata(:,2));
        permangleMat(isubj,ibin) = mean(permdata(:,1));
        permerrorMat(isubj,ibin) = mean(permdata(:,2));
    end
    
end

meandataangle = mean(angleMat);
meandataerror = mean(errorMat);
meanpermangle = mean(permangleMat);
meanpermerror = mean(permerrorMat);

semdataangle = std(angleMat)./sqrt(nSubj);
semdataerror = std(errorMat)./sqrt(nSubj);
sempermangle = std(permangleMat)./sqrt(nSubj);
sempermerror = std(permerrorMat)./sqrt(nSubj);

% plot(angleMat',errorMat'); hold on
% plot(mean(angleMat),mean(errorMat),'k-','LineWidth',3)

figure; hold on
fill([meandataangle fliplr(meandataangle)],[meandataerror-semdataerror fliplr(meandataerror+semdataerror)],[1 0.8 0.1],'FaceAlpha',0.5,'EdgeColor','none')
fill([meanpermangle fliplr(meanpermangle)],[meanpermerror-sempermerror fliplr(meanpermerror+sempermerror)],'k','FaceAlpha',0.2,'EdgeColor','none')
plot([90 90; 180 180; 270 270],[1.5 2.6],'k--')
ylim([1.5 2.6])

figure; hold on
fill([(meanpermerror+sempermerror).*cosd(meanpermangle) (meanpermerror(1)+sempermerror(1)).*cosd(meanpermangle(1)) ,...
    fliplr([(meanpermerror-sempermerror).*cosd(meanpermangle) (meanpermerror(1)-sempermerror(1)).*cosd(meanpermangle(1))])],...
    [(meanpermerror+sempermerror).*sind(meanpermangle) (meanpermerror(1)+sempermerror(1)).*sind(meanpermangle(1)),...
    fliplr([(meanpermerror-sempermerror).*sind(meanpermangle) (meanpermerror(1)-sempermerror(1)).*sind(meanpermangle(1))])],'k','FaceAlpha',0.2,'EdgeColor','none')
fill([(meandataerror+semdataerror).*cosd(meandataangle) (meandataerror(1)+semdataerror(1)).*cosd(meandataangle(1)) ,...
    fliplr([(meandataerror-semdataerror).*cosd(meandataangle) (meandataerror(1)-semdataerror(1)).*cosd(meandataangle(1))])],...
    [(meandataerror+semdataerror).*sind(meandataangle) (meandataerror(1)+semdataerror(1)).*sind(meandataangle(1)),...
    fliplr([(meandataerror-semdataerror).*sind(meandataangle) (meandataerror(1)-semdataerror(1)).*sind(meandataangle(1))])],[1 0.8 0.1],'FaceAlpha',0.5,'EdgeColor','none')
plot([-2.5 0; 2.5 0],[0 -2.5; 0 2.5],'k-')
axis equal
axis([-2.5 2.5 -2.5 2.5])
defaultplot

%% PERMUTATION TEST FOR STIMULUS: EXP 2 (SEPARATED PRIORITIES)
clear all
% close all

nBins = 20;
expnumber = 2; 

switch expnumber
    case 1
        load('exp1_zuzprocesseddata.mat')
%         group_data = group_data(:,[1 6 7 13]);
    case 2
        load('exp2_zuzprocesseddata.mat')
        group_data = group_data(:,[1 2 6 7 13]);
        subjVec = 4:14;
        dataa = sqrt((10.*cosd(group_data(:,5)) - group_data(:,3)).^2 +...
            (10.*sind(group_data(:,5)) - group_data(:,4)).^2);
        group_data = [group_data(:,[1 2 5]) dataa];
end
nSubj = length(subjVec);
priorityVec = [0.6 0.3 0.1];

group_data(any(isnan(group_data), 2), :) = [];

[angleMat, errorMat, permerrorMat, permangleMat] = deal(cell(1,3));

for isubj = 1:nSubj;
    subjnum = subjVec(isubj);
    
    for ipriority = 1:3;
        if isubj == 1;
        [angleMat{ipriority}, errorMat{ipriority}, permerrorMat{ipriority},...
            permangleMat{ipriority}] = deal(nan(nSubj,nBins));
        end
        
        priority = priorityVec(ipriority);
        
        idx = (group_data(:,1) == subjnum) & (group_data(:,2) == priority);
        nTrials = sum(idx);
        trialss = round(linspace(1,nTrials,nBins+1));
        
        dataa = group_data(idx,3:4);
        dataa = sortrows(dataa,1);
        
        permuteidx = randperm(nTrials);
        permutedataa = [dataa(:,1) dataa(permuteidx,2)];
        
        for ibin = 1:nBins;
            bindata = dataa(trialss(ibin):trialss(ibin+1),:);
            permdata = permutedataa(trialss(ibin):trialss(ibin+1),:);
            
            angleMat{ipriority}(isubj,ibin) = mean(bindata(:,1));
            errorMat{ipriority}(isubj,ibin) = mean(bindata(:,2));
            permangleMat{ipriority}(isubj,ibin) = mean(permdata(:,1));
            permerrorMat{ipriority}(isubj,ibin) = mean(permdata(:,2));
        end
    end
    
end

colorMat = {'r','b','k'};
figure; hold on
for ipriority = 1:3;
meandataangle = mean(angleMat{ipriority});
meandataerror = mean(errorMat{ipriority});
meanpermangle = mean(permangleMat{ipriority});
meanpermerror = mean(permerrorMat{ipriority});

semdataangle = std(angleMat{ipriority})./sqrt(nSubj);
semdataerror = std(errorMat{ipriority})./sqrt(nSubj);
sempermangle = std(permangleMat{ipriority})./sqrt(nSubj);
sempermerror = std(permerrorMat{ipriority})./sqrt(nSubj);

fill([(meanpermerror+sempermerror).*cosd(meanpermangle) (meanpermerror(1)+sempermerror(1)).*cosd(meanpermangle(1)) ,...
    fliplr([(meanpermerror-sempermerror).*cosd(meanpermangle) (meanpermerror(1)-sempermerror(1)).*cosd(meanpermangle(1))])],...
    [(meanpermerror+sempermerror).*sind(meanpermangle) (meanpermerror(1)+sempermerror(1)).*sind(meanpermangle(1)),...
    fliplr([(meanpermerror-sempermerror).*sind(meanpermangle) (meanpermerror(1)-sempermerror(1)).*sind(meanpermangle(1))])],'k','FaceAlpha',0.0,'EdgeColor',colorMat{ipriority})
fill([(meandataerror+semdataerror).*cosd(meandataangle) (meandataerror(1)+semdataerror(1)).*cosd(meandataangle(1)) ,...
    fliplr([(meandataerror-semdataerror).*cosd(meandataangle) (meandataerror(1)-semdataerror(1)).*cosd(meandataangle(1))])],...
    [(meandataerror+semdataerror).*sind(meandataangle) (meandataerror(1)+semdataerror(1)).*sind(meandataangle(1)),...
    fliplr([(meandataerror-semdataerror).*sind(meandataangle) (meandataerror(1)-semdataerror(1)).*sind(meandataangle(1))])],colorMat{ipriority},'FaceAlpha',0.3,'EdgeColor','none')

end
plot([-4 0; 4 0],[0 -4; 0 4],'k-')
axis equal
axis([-4 4 -4 4])



