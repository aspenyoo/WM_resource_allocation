
clear all

nPriorities = 3;
colorMat = [0.5 0.5 0.5; 0.5 0.5 1; 1 0.5 0.5 ];
markerMat = {'.','x','s'};

modelnameVec = {'optimal','free','fixed'}'

%% =================== EXP 1 ==========================
%           MAIN BEHAVIORAL RESULTS
% =====================================================

figure; hold on

expnumber = 1;
load('cleandata_nodisc.mat','data')
nSubj = 14;

% get data
filepath = ['fits/exp' num2str(expnumber) '/'];
dataerrorVec = nan(nSubj,nPriorities);
for ipriority = 1:nPriorities
    for isubj = 1:nSubj
        dataerrorVec(isubj,ipriority) = mean(data{isubj}{ipriority}(:,1));
    end
end

% M and SEM
dataMerror = mean(dataerrorVec);
dataSEMerror = std(dataerrorVec)./sqrt(nSubj);


for imodel = 2:3
    load([filepath 'modelpred_exp' num2str(expnumber) '_model' num2str(imodel) '.mat'],'preddata')
    
    modelerrorVec= nan(nSubj,nPriorities);
    for ipriority = 1:nPriorities
        for isubj = 1:nSubj
            modelerrorVec(isubj,ipriority) = mean(preddata{isubj}{ipriority}(:,1));
        end
    end
    
    % M and SEM
    modelMerror = mean(modelerrorVec);
    modelSEMerror = std(modelerrorVec)./sqrt(nSubj);
    
    
    dx = 0.2;
    %  ========== plot ==========
    subplot(1,2,imodel-1); hold on;
    for ipriority = 1:nPriorities
        mm = modelMerror(nPriorities+1-ipriority);
        sem = modelSEMerror(nPriorities+1-ipriority);
        fill([ipriority-dx ipriority+dx ipriority+dx ipriority-dx],...
            [mm-sem mm-sem mm+sem mm+sem],...
            colorMat(nPriorities+1-ipriority,:),'EdgeColor','none','FaceAlpha',0.4);
        errorbar(ipriority,dataMerror(nPriorities+1-ipriority),dataSEMerror(nPriorities+1-ipriority),...
            'k','Marker',markerMat{ipriority},'LineWidth',1);
    end
    defaultplot;
    axis([0.5 3.5 1 2.5])
    set(gca,'YTick',1:.5:2.5,'YTickLabel',{1, '','', 2.5});
    xlabel('priority'); ylabel('error');
    set(gca,'XTick',1:3,'XTickLabel',{'low','medium','high'})
    title([modelnameVec{imodel} ' model'])
    
end

%% % % % %  MODELING RESULTS % % % % % % % %





%% =================== EXP 2 ==========================
%           MAIN BEHAVIORAL RESULTS
% =====================================================

figure; hold on

expnumber = 2;
load('cleandata.mat','data')
nSubj = 11;

% load simulated data
filepath = ['fits/exp' num2str(expnumber) '/'];
[dataerrorVec, datadiscsizeVec ] = deal(nan(nSubj,nPriorities));
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

for imodel = 1:3
    load([filepath 'modelpred_exp' num2str(expnumber) '_model' num2str(imodel) '.mat'],'preddata')
    
    [modelerrorVec, modeldiscsizeVec ] = deal(nan(nSubj,nPriorities));
    for ipriority = 1:nPriorities
        for isubj = 1:nSubj
            modelerrorVec(isubj,ipriority) = mean(preddata{isubj}{ipriority}(:,1));
            modeldiscsizeVec(isubj,ipriority) = mean(preddata{isubj}{ipriority}(:,2));
        end
    end
    
    %  ========= calculate M, SEM across subjects ========
    modelMerror = mean(modelerrorVec);
    modelSEMerror = std(modelerrorVec)./sqrt(nSubj);
    
    modelMdiscsize = mean(modeldiscsizeVec);
    modelSEMdiscsize = std(modeldiscsizeVec)/sqrt(nSubj);
    
    
    subplot(3,2,2*imodel-1);hold on;
    for ipriority = 1:nPriorities
        mm = modelMerror(nPriorities+1-ipriority);
        sem = modelSEMerror(nPriorities+1-ipriority);
        fill([ipriority-dx ipriority+dx ipriority+dx ipriority-dx],...
            [mm-sem mm-sem mm+sem mm+sem],...
            colorMat(nPriorities+1-ipriority,:),'EdgeColor','none','FaceAlpha',0.4);
        errorbar(ipriority,dataMerror(nPriorities+1-ipriority),dataSEMerror(nPriorities+1-ipriority),...
            'k','Marker',markerMat{ipriority},'LineWidth',1);
    end
    defaultplot;
    axis([0.5 3.5 1 4])
    xlabel('priority'); ylabel('error');
    set(gca,'YTick',1:.4:2.2, 'YTickLabel', {1,'','', 2.2},...
        'XTick',1:3,'XTickLabel',{'low','medium','high'})
    
    
    subplot(3,2,2*imodel); hold on;
    for ipriority = 1:nPriorities
        mm = modelMdiscsize(nPriorities+1-ipriority);
        sem = modelSEMdiscsize(nPriorities+1-ipriority);
        fill([ipriority-dx ipriority+dx ipriority+dx ipriority-dx],...
            [mm-sem mm-sem mm+sem mm+sem],...
            colorMat(nPriorities+1-ipriority,:),'EdgeColor','none','FaceAlpha',0.4);
        errorbar(ipriority,dataMdiscsize(nPriorities+1-ipriority),dataSEMdiscsize(nPriorities+1-ipriority),...
            'k','Marker',markerMat{ipriority},'LineWidth',1);
    end
    defaultplot;
    axis([0.5 3.5 2 5])
    xlabel('priority'); ylabel('disc size')
    set(gca,'YTick',2:4,'YTickLabel',{2,'', 4},...
        'XTick',1:3,'XTickLabel',{'low','medium','high'})
end



%% % % % %    MODELING RESULTS % % % % % % % %

