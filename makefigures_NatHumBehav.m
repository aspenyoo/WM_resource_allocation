


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
colorMat = [0.5 0.5 0.5; 0.5 0.5 1; 1 0.5 0.5];
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
axis([0.5 3.5 1.6 2.4])
set(gca,'YTick',1.6:.4:2.4,'YTickLabel',1.6:.4:2.4);

for ipriority = 1:nPriorities
    errorbar(ipriority,dataMerror(nPriorities+1-ipriority),dataSEMerror(nPriorities+1-ipriority),...
        'k','Marker',markerMat{ipriority},'LineWidth',1);
    errorb(ipriority,dataMerror(nPriorities+1-ipriority),dataSEMerror(nPriorities+1-ipriority),...
        'LineWidth',1);
end


subplot(1,2,2); hold on;
defaultplot;
axis([0.5 3.5 2.2 3.8])
set(gca,'YTick',[2.2 3 3.8],'YTickLabel',[2.2 3 3.8]);
xlabel('priority'); ylabel('error');
set(gca,'XTick',1:3,'XTickLabel',{'low','medium','high'})
title('disc size')
for ipriority = 1:nPriorities
    errorbar(ipriority,dataMdiscsize(nPriorities+1-ipriority),dataSEMdiscsize(nPriorities+1-ipriority),...
        'k','Marker',markerMat{ipriority},'LineWidth',1);
    errorb(ipriority,dataMdiscsize(nPriorities+1-ipriority),dataSEMdiscsize(nPriorities+1-ipriority),...
        'LineWidth',1);
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
colorMat = [0.3 0.3 0.3; 0.3 0.3 1; 1 0.3 0.3];

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
set(hter,'marker','o','markerfacecolor','k','markersize',4','markeredgecolor','k')

%% ====================================================
%          MODELING RESULTS: COMBINED MODELS
% =====================================================

clear all
expnumber = 2;
fixedrisk = [];%'_fixedrisk';
% colorMat = [aspencolors('dustyrose'); aspencolors('booger'); aspencolors('seacolored')];
colorMat = [0 0 1; 1 0 0; 0 0 0];

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
    ha = tight_subplot(3,3,{[.03 .03],[.03 .07]},[.1 .01],[.1 .01]);
    set(gcf,'Position',[28 504 800 500])
else
    ha = tight_subplot(1,3,.03,[.26 .05],[.11 .05]);
    set(gcf,'Position',[28 504 800 236])
end

% load and plot predicted data for each model
for imodel = (3-expnumber):3
    
    % load data
    load([filepath 'modelpred_exp' num2str(expnumber) '_model' num2str(imodel) fixedrisk '.mat'],'preddata')
    
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
    
    for ipriority = 1:nPriorities
        
        if(expnumber == 2)
            axes(ha(3*ipriority-2))
        else
            axes(ha(ipriority))
        end
        fill([xlims fliplr(xlims)],[meansimerror{ipriority}-semsimerror{ipriority}...
            fliplr(meansimerror{ipriority}+semsimerror{ipriority})],colorMat(imodel,:),'EdgeColor','none','FaceAlpha',0.4);
        hold on;
        
        if (expnumber == 2)
            % discsize
            axes(ha(3*ipriority-1))
            fill([xlims fliplr(xlims)],[meansimdiscsize{ipriority}-semsimdiscsize{ipriority}...
                fliplr(meansimdiscsize{ipriority}+semsimdiscsize{ipriority})],colorMat(imodel,:),'EdgeColor','none','FaceAlpha',0.4);
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
        
        % ================ group plot ====================
        
        meanmeanquantsimerror = cellfun(@nanmean,meanquantsimerror,'UniformOutput',false);
        meanmeanquantsimdiscsize = cellfun(@mean,meanquantsimdiscsize,'UniformOutput',false);
        semmeanquantsimdiscsize= cellfun(@(x) std(x)./sqrt(size(x,1)),meanquantsimdiscsize,'UniformOutput',false);
        
        for ipriority = 1:nPriorities
            axes(ha(3*ipriority))
            %         subplot(3,3,6+ipriority)
            hold on
            plot_summaryfit(meanmeanquantsimerror{ipriority},[],[],meanmeanquantsimdiscsize{ipriority},...
                semmeanquantsimdiscsize{ipriority},[],colorMat(imodel,:))
        end
    end
    
end

% ============= PLOT REAL DATA ======================

% load data
load(['exp' num2str(expnumber) '_cleandata.mat'],'data')
load([filepath 'fits_model' num2str(imodel) '.mat'])


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

% figure;
% colorMat = {'r','b','k'};
% if (expnumber == 2)
%     ha = tight_subplot(3,3,{[.03 .03],[.03 .07]},[.1 .01],[.1 .01]);
%     set(gcf,'Position',[28 504 800 500])
% else
%     ha = tight_subplot(1,3,.03,[.26 .05],[.11 .05]);
%     set(gcf,'Position',[28 504 800 236])
% end


meanerror = cellfun(@mean,error,'UniformOutput',false);
semerror = cellfun(@(x) std(x)./sqrt(size(x,1)),error,'UniformOutput',false);

if (expnumber == 2)
    meandiscsize = cellfun(@mean,discsize,'UniformOutput',false);
    semdiscsize = cellfun(@(x) std(x)./sqrt(size(x,1)),discsize,'UniformOutput',false);
end

for ipriority = 1:nPriorities
    
    if(expnumber == 2)
        axes(ha(3*ipriority-2))
    else
        axes(ha(ipriority))
    end
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
    
    % ================ group plot ====================
    
    meanmeanquanterror = cellfun(@mean,meanquanterror,'UniformOutput',false);
    meanmeanquantdiscsize = cellfun(@mean,meanquantdiscsize,'UniformOutput',false);
    semmeanquantdiscsize= cellfun(@(x) std(x)./sqrt(size(x,1)),meanquantdiscsize,'UniformOutput',false);

    for ipriority = 1:nPriorities
        axes(ha(3*ipriority))
        %         subplot(3,3,6+ipriority)
        hold on
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

%% ====================================================
%          MODELING RESULTS: COMBINED PRIORITIES
% =====================================================

clear all
expnumber = 2;
fixedrisk = [];%'_fixedrisk';
% colorMat = [aspencolors('dustyrose'); aspencolors('booger'); aspencolors('seacolored')];
colorMat = [0 0 1; 1 0 0; 0 0 0];

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
    ha = tight_subplot(3,3,{[.03 .03],[.03 .07]},[.1 .01],[.1 .01]);
    set(gcf,'Position',[28 504 800 500])
else
    ha = tight_subplot(1,2,.03,[.05],[.05]);
    set(gcf,'Position',[28 504 800 236])
end

% load and plot predicted data for each model
for imodel = (3-expnumber):3
    
    % load data
    load([filepath 'modelpred_exp' num2str(expnumber) '_model' num2str(imodel) fixedrisk '.mat'],'preddata')
    
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
    
    for ipriority = 1:nPriorities
        
        if(expnumber == 2)
            axes(ha(3*imodel-2))
        else
            axes(ha(imodel-1))
        end
        fill([xlims fliplr(xlims)],[meansimerror{ipriority}-semsimerror{ipriority}...
            fliplr(meansimerror{ipriority}+semsimerror{ipriority})],colorMat(ipriority,:),'EdgeColor','none','FaceAlpha',0.4);
        hold on;
        
        if (expnumber == 2)
            % discsize
            axes(ha(3*imodel-1))
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
        
        % ================ group plot ====================
        
        meanmeanquantsimerror = cellfun(@nanmean,meanquantsimerror,'UniformOutput',false);
        meanmeanquantsimdiscsize = cellfun(@mean,meanquantsimdiscsize,'UniformOutput',false);
        semmeanquantsimdiscsize= cellfun(@(x) std(x)./sqrt(size(x,1)),meanquantsimdiscsize,'UniformOutput',false);
        
        for ipriority = 1:nPriorities
            axes(ha(3*imodel))
            %         subplot(3,3,6+ipriority)
            hold on
            plot_summaryfit(meanmeanquantsimerror{ipriority},[],[],meanmeanquantsimdiscsize{ipriority},...
                semmeanquantsimdiscsize{ipriority},[],colorMat(ipriority,:))
        end
    end
    
end

% ============= PLOT REAL DATA ======================

% load data
load(['exp' num2str(expnumber) '_cleandata.mat'],'data')
load([filepath 'fits_model' num2str(imodel) '.mat'])


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

% figure;
% colorMat = {'r','b','k'};
% if (expnumber == 2)
%     ha = tight_subplot(3,3,{[.03 .03],[.03 .07]},[.1 .01],[.1 .01]);
%     set(gcf,'Position',[28 504 800 500])
% else
%     ha = tight_subplot(1,3,.03,[.26 .05],[.11 .05]);
%     set(gcf,'Position',[28 504 800 236])
% end
for imodel = (3-expnumber):3

meanerror = cellfun(@mean,error,'UniformOutput',false);
semerror = cellfun(@(x) std(x)./sqrt(size(x,1)),error,'UniformOutput',false);

if (expnumber == 2)
    meandiscsize = cellfun(@mean,discsize,'UniformOutput',false);
    semdiscsize = cellfun(@(x) std(x)./sqrt(size(x,1)),discsize,'UniformOutput',false);
end

for ipriority = 1:nPriorities
    
    if(expnumber == 2)
        axes(ha(3*imodel-2))
    else
        axes(ha(imodel-1))
    end
    hold on;
    errorbar(xlims,meanerror{ipriority},semerror{ipriority},'Color',colorMat(ipriority,:),'LineStyle','none','LineWidth',1);
    defaultplot
    axis([0 10 0 0.4])
    
    if expnumber == 2
        %          axis([0 10 0 0.6])
        if imodel == 3
            xlabel('error');
%         else
%             set(ha(3*imodel-2),'XTickLabel','');
        end
        ylabel('proportion','FontSize',14);
        
    else
        
        xlabel('error','FontSize',16)
        set(ha(imodel-1),'YTick',[0 0.2 0.4],'FontSize',12);
        
        if imodel == 2
            ylabel('proportion','FontSize',16)
            
        else
            set(ha(imodel-1),'YTickLabel','');
            
        end
    end
    
    
    if (expnumber == 2)
        % discsize
        axes(ha(3*imodel-1))
        hold on;
        errorbar(xlims,meandiscsize{ipriority},semdiscsize{ipriority},'Color',colorMat(ipriority,:),'LineStyle','none','LineWidth',1);
        defaultplot
        axis([0 10 0 0.6])
        if imodel == 3
            xlabel('disc size','FontSize',14);
%         else
%             set(ha(3*imodel-1),'XTickLabel','');
        end
        set(ha(3*imodel-1),'YTickLabel','');
        
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
    
    % ================ group plot ====================
    
    meanmeanquanterror = cellfun(@mean,meanquanterror,'UniformOutput',false);
    meanmeanquantdiscsize = cellfun(@mean,meanquantdiscsize,'UniformOutput',false);
    semmeanquantdiscsize= cellfun(@(x) std(x)./sqrt(size(x,1)),meanquantdiscsize,'UniformOutput',false);

    for ipriority = 1:nPriorities
        axes(ha(3*imodel))
        %         subplot(3,3,6+ipriority)
        hold on
        plot_summaryfit(meanmeanquanterror{ipriority},meanmeanquantdiscsize{ipriority},semmeanquantdiscsize{ipriority},...
            [],[],colorMat(ipriority,:))
        
        ylabel('disc size');
        axis([0 6 2 6])
        set(ha(3*imodel),'XTick',[0 3 6],'YTick',[2 4 6]);
        set(ha(3*imodel),'YTickLabel',[2 4 6])
        if (expnumber == 1)
            set(ha(imodel),'XTickLabel',[0 3 6])
            xlabel('error');
        end
    end
    set(ha(3*imodel),'XTickLabel',[0 3 6]);
    xlabel('error');
end
end