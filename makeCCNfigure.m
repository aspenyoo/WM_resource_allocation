clear all

imodel = 2;
fixedrisk = [];%'_fixedrisk';
nPriorities = 3;
% colorMat = [ .5 .5 1; 1 1 .5; 1 .5 .5];
% colorMat = [aspencolors('indigo'); aspencolors('bluegreen'); aspencolors('orange')];
colorMat = [0.5 0.5 0.5; 0.5 0.5 1; 1 0.5 0.5 ];
markerMat = {'.','x','s'};

figure;
set(gcf,'Position',[28 504 1300 236])

% ========================
% EXP 1
% ========================
expnumber = 1;
load('cleandata_nodisc.mat','data')
nSubj = 14;

filepath = ['fits/exp' num2str(expnumber) fixedrisk '/'];
load([filepath 'modelpred_exp' num2str(expnumber) '_model' num2str(imodel) fixedrisk '.mat'],'preddata')

[dataerrorVec, modelerrorVec ] = deal(nan(nSubj,nPriorities));
for ipriority = 1:nPriorities
    
    for isubj = 1:nSubj
        dataerrorVec(isubj,ipriority) = mean(preddata{isubj}{ipriority}(:,1));
        modelerrorVec(isubj,ipriority) = mean(preddata{isubj}{ipriority}(:,1));
    end
    
end

%  ========= calculate M, SEM across subjects ========
dataMerror = mean(dataerrorVec);
dataSEMerror = std(dataerrorVec)./sqrt(nSubj);

modelMerror = mean(modelerrorVec);
modelSEMerror = std(modelerrorVec)./sqrt(nSubj);


dx = 0.2;
%  ========== plot ==========
% subplot(1,9,1); hold on;
% for ipriority = 1:nPriorities
%     mm = modelMerror(nPriorities+1-ipriority);
%     sem = modelSEMerror(nPriorities+1-ipriority);
%     fill([ipriority-dx ipriority+dx ipriority+dx ipriority-dx],...
%         [mm-sem mm-sem mm+sem mm+sem],...
%         colorMat(nPriorities+1-ipriority,:),'EdgeColor','none','FaceAlpha',0.4);
%     errorbar(ipriority,dataMerror(nPriorities+1-ipriority),dataSEMerror(nPriorities+1-ipriority),...
%         'k','Marker',markerMat{ipriority},'LineWidth',1); 
% end
% defaultplot;
% axis([0.5 3.5 1 2.2])
% set(gca,'YTick',1:.4:2.2,'YTickLabel',{1, '','', 2.2});
% xlabel('priority'); ylabel('error');
% set(gca,'XTick',1:3,'XTickLabel',{'low','medium','high'})

% plot bargraph instead of errorbars
dx = 0.1;
subplot(1,3,1);hold on;
axis([0.5 3.5 1 2.2])
set(gca,'YTick',1:.4:2.2,'YTickLabel',{1, '','', 2.2});
xlabel('priority'); ylabel('error');
set(gca,'XTick',1:3,'XTickLabel',{'low','medium','high'})
for ipriority = 1:nPriorities
    mm = modelMerror(nPriorities+1-ipriority);
    sem = modelSEMerror(nPriorities+1-ipriority);
    errorbar(ipriority-dx,dataMerror(nPriorities+1-ipriority),dataSEMerror(nPriorities+1-ipriority),...
        'k','Marker',markerMat{ipriority},'LineWidth',1); 
    errorbar(ipriority+dx,modelMerror(nPriorities+1-ipriority),modelSEMerror(nPriorities+1-ipriority),...
        'Color',colorMat(ipriority,:),'LineWidth',1); 
end
defaultplot;


% ------------ ternary plot ------------
load([filepath 'fits_model' num2str(imodel) '.mat'])

pMat = ML_parameters(:,end-1:end);
pMat(:,3) = 1-sum(pMat,2);

subplot(1,9,6:7);
% axis
[h,hg,htick]=terplot;
c1 = [1 0.5 1/3];
c2 = [0 0.5 1/3];
c3 = [0 0 1/3];
hlabels=terlabel('high','medium','low');

% colors
axis1 = colorMat(2,:);
axis2 = colorMat(1,:);
axis3 = colorMat(3,:);

%
set(h,'LineWidth',1)

%--  Change the color of the grid lines
set(hg(:,1),'color',axis1)
set(hg(:,2),'color',axis2)
set(hg(:,3),'color',axis3)

% make 0.6 0.3 0.1 solid
set(hg(3,1),'LineStyle','--')
set(hg(6,2),'LineStyle','--')
set(hg(1,3),'LineStyle','--')

%--  Modify the labels
set(hlabels,'fontsize',12)
set(hlabels(1),'color',axis1)
set(hlabels(2),'color',axis2)
set(hlabels(3),'color',axis3)

%--  Modify the tick labels
set(htick(:,1),'color',axis1,'linewidth',3)
set(htick(:,2),'color',axis2,'linewidth',3)
set(htick(:,3),'color',axis3,'linewidth',3)

% plot data
hter=ternaryc(pMat(:,1),pMat(:,2),pMat(:,3));
set(hter,'marker','o','markerfacecolor','k','markersize',4','markeredgecolor','k')




% =======================
% EXP 2
% =======================
expnumber = 2;
load('cleandata.mat','data')
nSubj = 11;


% load simulated data
filepath = ['fits/exp' num2str(expnumber) fixedrisk '/'];
load([filepath 'modelpred_exp' num2str(expnumber) '_model' num2str(imodel) fixedrisk '.mat'],'preddata')

[dataerrorVec, datadiscsizeVec, modelerrorVec, modeldiscsizeVec ] = deal(nan(nSubj,nPriorities));
for ipriority = 1:nPriorities
    
    for isubj = 1:nSubj
        dataerrorVec(isubj,ipriority) = mean(preddata{isubj}{ipriority}(:,1));
        datadiscsizeVec(isubj,ipriority) = mean(preddata{isubj}{ipriority}(:,2));
        modelerrorVec(isubj,ipriority) = mean(preddata{isubj}{ipriority}(:,1));
        modeldiscsizeVec(isubj,ipriority) = mean(preddata{isubj}{ipriority}(:,2));
    end
    
end

%  ========= calculate M, SEM across subjects ========
dataMerror = mean(dataerrorVec);
dataSEMerror = std(dataerrorVec)./sqrt(nSubj);

dataMdiscsize = mean(datadiscsizeVec);
dataSEMdiscsize = std(datadiscsizeVec)/sqrt(nSubj);

modelMerror = mean(modelerrorVec);
modelSEMerror = std(modelerrorVec)./sqrt(nSubj);

modelMdiscsize = mean(modeldiscsizeVec);
modelSEMdiscsize = std(modeldiscsizeVec)/sqrt(nSubj);

% plot modelpredictions as errorbars

%  ========== plot ==========
% subplot(1,9,2);hold on;
% for ipriority = 1:nPriorities
%     mm = modelMerror(nPriorities+1-ipriority);
%     sem = modelSEMerror(nPriorities+1-ipriority);
%     fill([ipriority-dx ipriority+dx ipriority+dx ipriority-dx],...
%         [mm-sem mm-sem mm+sem mm+sem],...
%         colorMat(nPriorities+1-ipriority,:),'EdgeColor','none','FaceAlpha',0.4);
%     errorbar(ipriority,dataMerror(nPriorities+1-ipriority),dataSEMerror(nPriorities+1-ipriority),...
%         'k','Marker',markerMat{ipriority},'LineWidth',1); 
% end
% defaultplot;
% axis([0.5 3.5 1 2.2])
% xlabel('priority'); ylabel('error');
% set(gca,'YTick',1:.4:2.2, 'YTickLabel', {1,'','', 2.2},...
%     'XTick',1:3,'XTickLabel',{'low','medium','high'})
% 
% 
% subplot(1,9,3); hold on;
% for ipriority = 1:nPriorities
%     mm = modelMdiscsize(nPriorities+1-ipriority);
%     sem = modelSEMdiscsize(nPriorities+1-ipriority);
%     fill([ipriority-dx ipriority+dx ipriority+dx ipriority-dx],...
%         [mm-sem mm-sem mm+sem mm+sem],...
%         colorMat(nPriorities+1-ipriority,:),'EdgeColor','none','FaceAlpha',0.4);
%     errorbar(ipriority,dataMdiscsize(nPriorities+1-ipriority),dataSEMdiscsize(nPriorities+1-ipriority),...
%         'k','Marker',markerMat{ipriority},'LineWidth',1); 
% end
% defaultplot;
% axis([0.5 3.5 2 4])
% xlabel('priority'); ylabel('disc size')
% set(gca,'YTick',2:4,'YTickLabel',{2,'', 4},...
%     'XTick',1:3,'XTickLabel',{'low','medium','high'})

subplot(1,3,2);hold on;
for ipriority = 1:nPriorities
%     fill([ipriority-dx ipriority+dx ipriority+dx ipriority-dx],...
%         [mm-sem mm-sem mm+sem mm+sem],...
%         colorMat(nPriorities+1-ipriority,:),'EdgeColor','none','FaceAlpha',0.4);
    errorb(ipriority-dx,dataMerror(nPriorities+1-ipriority),dataSEMerror(nPriorities+1-ipriority),...
        'k','Marker',markerMat{ipriority},'LineWidth',1); 
    errorb(ipriority+dx,modelMerror(nPriorities+1-ipriority),modelSEMerror(nPriorities+1-ipriority),...
        'Color',colorMat(ipriority,:),'LineWidth',1); 
end
defaultplot;
axis([0.5 3.5 1 2.2])
xlabel('priority'); ylabel('error');
set(gca,'YTick',1:.4:2.2, 'YTickLabel', {1,'','', 2.2},...
    'XTick',1:3,'XTickLabel',{'low','medium','high'})


subplot(1,3,3); hold on;
for ipriority = 1:nPriorities
    errorb(ipriority-dx,dataMdiscsize(nPriorities+1-ipriority),dataSEMdiscsize(nPriorities+1-ipriority),...
        'k','Marker',markerMat{ipriority},'LineWidth',1); 
    errorb(ipriority+dx,modelMdiscsize(nPriorities+1-ipriority),modelSEMdiscsize(nPriorities+1-ipriority),...
        'Color',colorMat(ipriority,:),'LineWidth',1); 
end
defaultplot;
axis([0.5 3.5 2 4])
xlabel('priority'); ylabel('disc size')
set(gca,'YTick',2:4,'YTickLabel',{2,'', 4},...
    'XTick',1:3,'XTickLabel',{'low','medium','high'})

% ========== quantile correlation plot per subject ===========

nQuants = 6;
for isubj = 1:11
    
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

    end
    
end

% ================ group plot ====================

meanmeanquanterror = cellfun(@mean,meanquanterror,'UniformOutput',false);
meanmeanquantdiscsize = cellfun(@mean,meanquantdiscsize,'UniformOutput',false);
semmeanquantdiscsize= cellfun(@(x) std(x)./sqrt(size(x,1)),meanquantdiscsize,'UniformOutput',false);
meanmeanquantsimerror = cellfun(@nanmean,meanquantsimerror,'UniformOutput',false);
meanmeanquantsimdiscsize = cellfun(@mean,meanquantsimdiscsize,'UniformOutput',false);
semmeanquantsimdiscsize= cellfun(@(x) std(x)./sqrt(size(x,1)),meanquantsimdiscsize,'UniformOutput',false);

for ipriority = 1:nPriorities
    subplot(1,9,4:5); hold on
    plot_summaryfit(meanmeanquantsimerror{ipriority},[],[],meanmeanquantsimdiscsize{ipriority},...
        semmeanquantsimdiscsize{ipriority},[],colorMat(ipriority,:));
end
for ipriority = 1:nPriorities
    subplot(1,9,4:5); hold on
    errorbar(meanmeanquanterror{nPriorities+1-ipriority},meanmeanquantdiscsize{nPriorities+1-ipriority},...
        semmeanquantsimdiscsize{nPriorities+1-ipriority},'k','Marker',markerMat{ipriority},...
        'LineStyle','none','LineWidth',1); 
end


    ylabel('disc size');
    axis([0 6 2 4])
    set(gca,'XTick',[0 3 6],'YTick',[2 3 4]);
    set(gca,'YTickLabel',[2 3 4])
    xlabel('error');
set(gca,'XTickLabel',[0 3 6]); 
xlabel('error');

% ------------ ternary plot ------------
load([filepath 'fits_model' num2str(imodel) '.mat'])

pMat = ML_parameters(:,end-1:end);
pMat(:,3) = 1-sum(pMat,2);

subplot(1,9,8:9);
% axis
[h,hg,htick]=terplot;
c1 = [1 0.5 1/3];
c2 = [0 0.5 1/3];
c3 = [0 0 1/3];
hlabels=terlabel('high','medium','low');

% colors
axis1 = colorMat(2,:);
axis2 = colorMat(1,:);
axis3 = colorMat(3,:);

%
set(h,'LineWidth',1)

%--  Change the color of the grid lines
set(hg(:,1),'color',axis1)
set(hg(:,2),'color',axis2)
set(hg(:,3),'color',axis3)

% make 0.6 0.3 0.1 solid
set(hg(3,1),'LineStyle','--')
set(hg(6,2),'LineStyle','--')
set(hg(1,3),'LineStyle','--')

%--  Modify the labels
set(hlabels,'fontsize',12)
set(hlabels(1),'color',axis1)
set(hlabels(2),'color',axis2)
set(hlabels(3),'color',axis3)

%--  Modify the tick labels
set(htick(:,1),'color',axis1,'linewidth',3)
set(htick(:,2),'color',axis2,'linewidth',3)
set(htick(:,3),'color',axis3,'linewidth',3)

% plot data
hter=ternaryc(pMat(:,1),pMat(:,2),pMat(:,3));
set(hter,'marker','o','markerfacecolor','k','markersize',4','markeredgecolor','k')
