

%% OBLIQUE EFFECT IN DATA?
% checking if there is any obvious stimulus-dependent noise in data

clear all; close all

expnumber = 2; 
collapseQuadrants = 1; 
makeDistanceFromCardinal = 0; 

switch expnumber
    case 1
        load('exp1_data_trialinfo.mat')
        subjVec = [1:9 11:15];
        
        finalerror = group_data(:,8);
        anglee = atand(group_data(:,14)./group_data(:,13));

        if (collapseQuadrants)
            anglee = abs(anglee);
        else
            idx = (group_data(:,13) < 0) & (group_data(:,14) > 0); % second quadrant
            anglee(idx) = anglee(idx) + 180;

            idx = (group_data(:,13) < 0) & (group_data(:,14) < 0); % third quadrant
            anglee(idx) = anglee(idx) + 180;

            idx = (group_data(:,13) > 0) & (group_data(:,14) < 0); % fourth quadrant
            anglee(idx) = anglee(idx) + 360;
        end

        group_data = [group_data(:,1) anglee finalerror]; % subject number, angle location, final error
    case 2
        load('exp2_data_trialinfo.mat')
        group_data = group_data(:,[1 6 7 13]);
        subjVec = 4:14;
        
        anglee = group_data(:,4);
        if (collapseQuadrants)
            idx = (anglee >= 90) & (anglee < 180); % second quadrant
            anglee(idx) = 180 - anglee(idx); 
            
            idx = (anglee >= 180) & (anglee < 270); % third quadrant
            anglee(idx) = anglee(idx) - 180; 
            
            idx = (anglee >= 270) & (anglee < 360); % fourth quadrant
            anglee(idx) = 360 - anglee(idx);
        end
        
        finalerror = sqrt((10.*cosd(group_data(:,4)) - group_data(:,2)).^2 +...
            (10.*sind(group_data(:,4)) - group_data(:,3)).^2);
        group_data = [group_data(:,1) anglee finalerror]; % subject number, angle location, final error
end
group_data(any(isnan(group_data), 2), :) = [];
if (makeDistanceFromCardinal)
    idx = group_data(:,2) > 45;
    group_data(idx,2) = 90 - group_data(idx,2);
end
nSubj = length(subjVec);

% correlation for each subject
[bVec, pVec] = deal(nan(1,nSubj));
bintVec = nan(nSubj,2);

angleVec = unique(group_data(:,2));
nAngles = length(angleVec);
mean_errors = nan(nSubj,nAngles);
for isubj = 1:nSubj;
    subjnum = subjVec(isubj);
    
    idx = group_data(:,1) == subjnum;
    nTrials = sum(idx);
    
    dataa = group_data(idx,2:3);
    [b,bint,~,~,stats] = regress(dataa(:,1), [ones(nTrials,1) dataa(:,2)]);
    bVec(isubj) = b(2);
    bintVec(isubj,:) = bint(2,:);
    pVec(isubj) = stats(3); 
    for iangle = 1:nAngles;
        angle = angleVec(iangle);
        idxx = dataa(:,1) == angle;
        
        mean_errors(isubj,iangle) = mean(dataa(idxx,2));
    end
end

figure; hold on
errorbar(angleVec,mean(mean_errors),std(mean_errors)./sqrt(nSubj),'k','LineStyle','none');
plot(angleVec,mean(mean_errors),'ko')
ylim([0 2.5])
defaultplot
%%


%% CORRELATION BETWEEN ERROR AND CIRCLE SIZE

clear all;

% load and preprocess data
% ------------------------
subjVec = 4:14;
priorityVec = [0.6 0.3 0.1];

load('exp2_data_trialinfo.mat')
finalerror = sqrt((10.*cosd(group_data(:,13)) - group_data(:,6)).^2 +...
    (10.*sind(group_data(:,13)) - group_data(:,7)).^2);

% collapsing all data onto first quadrant
collapsedangle = group_data(:,13);
idx = (collapsedangle < 180) & (collapsedangle > 90); % second quadrant
collapsedangle(idx) = 180 - collapsedangle(idx);
idx = (collapsedangle < 270) & (collapsedangle > 180); % third quadrant
collapsedangle(idx) = collapsedangle(idx) - 180;
idx = (collapsedangle < 360) & (collapsedangle > 270); % second quadrant
collapsedangle(idx) = 360 - collapsedangle(idx);

% columns correspond to [subject, priority, angle, error, circle size]
group_data = [group_data(:,1:2) collapsedangle finalerror group_data(:,4)];
group_data(any(isnan(group_data), 2), :) = [];
        
nSubj = length(subjVec);
nPriorities = length(priorityVec);

% put each subject's data in terms of z-scores
groupzscoredata = group_data;
[rMat, pMat, nTrials,Mat] = deal(nan(nSubj,nPriorities));
for isubj = 1:nSubj;
    subjnum = subjVec(isubj);
    
    idx = group_data(:,1) == subjnum;
    subjdata = group_data(idx,:);
    
    groupzscoredata(idx,4:5) = [zscore(subjdata(:,4)) zscore(subjdata(:,5))];

    for ipriority = 1:nPriorities;
        priority = priorityVec(ipriority);
        
        idxx = idx & (groupzscoredata(:,2) == priority);
        nTrialsMat(isubj,ipriority) = sum(idxx);
        [r,p] = corr(groupzscoredata(idxx,4:5),'Type','Spearman');
        rMat(isubj,ipriority) = r(2);
        pMat(isubj,ipriority) = p(2);
    end
end

[rVec, pVec] = deal(nan(1,nPriorities));
for ipriority = 1:nPriorities;
    priority = priorityVec(ipriority);
    idx = groupzscoredata(:,2) == priority;
    
    [r,p] = corr(groupzscoredata(idx,4:5),'Type','Spearman');
    rVec(ipriority) = r(2);
    pVec(ipriority) = p(2);
end

rMat
pMat
rVec
pVec
%% PERMUTATION TEST: STIMULUS LOCATION

clear all; clc

rng(1) % random seed for reproduceability

% load and preprocess data
% ------------------------
subjVec = 4:14;
priorityVec = [0.6 0.3 0.1];

load('exp2_data_trialinfo.mat')
finalerror = sqrt((10.*cosd(group_data(:,13)) - group_data(:,6)).^2 +...
    (10.*sind(group_data(:,13)) - group_data(:,7)).^2);

% collapsing all data onto first quadrant
collapsedangle = group_data(:,13);
idx = (collapsedangle < 180) & (collapsedangle > 90); % second quadrant
collapsedangle(idx) = 180 - collapsedangle(idx);
idx = (collapsedangle < 270) & (collapsedangle > 180); % third quadrant
collapsedangle(idx) = collapsedangle(idx) - 180;
idx = (collapsedangle < 360) & (collapsedangle > 270); % second quadrant
collapsedangle(idx) = 360 - collapsedangle(idx);

% columns correspond to [subject, priority, angle, error, circle size]
group_data = [group_data(:,1:2) collapsedangle finalerror group_data(:,4)];
group_data(any(isnan(group_data), 2), :) = [];
        
angleVec = unique(group_data(:,3));
nPerms = 1000;
nSubj = length(subjVec);
nAngles = length(angleVec);
nPriorities = length(priorityVec);

% put each subject's data in terms of z-scores
groupzscoredata = group_data;
for isubj = 1:nSubj;
    subjnum = subjVec(isubj);
    
    idx = group_data(:,1) == subjnum;
    subjdata = group_data(idx,:);
    
    groupzscoredata(idx,4:5) = [zscore(subjdata(:,4)) zscore(subjdata(:,5))];
end

% permute cicle size data for each bin for each subject
nTrialsVec = [];
permdata = nan(size(groupzscoredata,1),nPerms);
for isubj = 1:nSubj; % for each subject
    subjnum = subjVec(isubj);
    idx = groupzscoredata(:,1) == subjnum;
    
    for ipriority = 1:nPriorities % for each priority
        priority = priorityVec(ipriority);
        idxx = idx & (groupzscoredata(:,2) == priority);

        for iangle = 1:nAngles; % for each angle
            angle = angleVec(iangle);
            idxxx = idxx & (groupzscoredata(:,3) == angle); % indices of all trials that fit subject, priority, and angle criteria
            
            nTrials = sum(idxxx); % number of trials that fit criteria
            nTrialsVec = [nTrialsVec nTrials];
            % make a permutation matrix
            permMat = nan(nTrials,nPerms); % permuting these trials nPerms times
            for iperm = 1:nPerms; % for each permutation
                blah = randperm(nTrials);
                while any(blah == 1:nTrials) % implement derangement
                    blah = randperm(nTrials);
                end
                permMat(:,iperm) = blah';
            end
            
            circlesizedata = group_data(idxxx,5); % all wager data for the given subject, priority, and angle
            permdata(idxxx,:) = circlesizedata(permMat); % permuted data 1000 times
        end
    end
    
end

% calculate permuted null distributions
% -------------------------------------

[corrVec, pVec] = deal(nan(1,4*nSubj)); % real correlations and significancee values
[corrMat, pMat] = deal(nan(nPerms,4*nSubj)); % permuted values
for isubj = 1:nSubj; % for each subject
    subjnum = subjVec(isubj);
    idx = groupzscoredata(:,1) == subjnum;

    for ipriority = 1:nPriorities
        priority = priorityVec(ipriority);
         idxx = idx & (groupzscoredata(:,2) == priority);
         
         [r, p] = corrcoef([groupzscoredata(idxx,4:5), permdata(idxx,:)]);
         corrVec(4*(isubj-1) + ipriority) = r(2,1);
         pVec(4*(isubj-1) +ipriority) = p(2,1);
         corrMat(:,4*(isubj-1) +ipriority) = r(3:end,1);
         pMat(:,4*(isubj-1) +ipriority) = p(3:end,1);
    end
end

% print M \pm SEM real and permuated correlations
% -----------------------------------------------

% actual values
corrVec(4:4:44) = [];
data_M = mean(corrVec)
data_SEM = std(corrVec)/sqrt(length(corrVec))

% median of correlation distribution ofpermutations
meds = quantile(corrMat,.5);
meds(4:4:44) = [];
perm_M = mean(meds)
perm_SEM = std(meds)/sqrt(length(meds))


% wilcoxon signed-rank test
% -------------------------
[p,h,stats] = signrank(meds,corrVec)

%% MEAN EFFECT OF DELAY ON ERROR

clear all; close all; clc

% load and preprocess data
% ------------------------
subjVec = 4:14;
priorityVec = [0.6 0.3 0.1];

load('exp2_data_trialinfo.mat')
finalerror = sqrt((10.*cosd(group_data(:,13)) - group_data(:,6)).^2 +...
    (10.*sind(group_data(:,13)) - group_data(:,7)).^2);

% columns correspond to [subject, priority, delay, error, circle size]
group_data = [group_data(:,1:3) finalerror group_data(:,4)];
group_data(any(isnan(group_data), 2), :) = [];
   
delayVec = 1000:500:4000;
nSubj = length(subjVec);
nDelays = length(delayVec);
nPriorities = length(priorityVec);

% combine dleays that are very close to one another
for idelay = 1:nDelays;
    delay = delayVec(idelay);
    
    idx = (delay-100 < group_data(:,3)) & (delay+100 > group_data(:,3));
    group_data(idx,3) = delay;
end

[bVec, pVec] = deal(nan(1,nSubj));
bintVec = nan(nSubj,2);
errorMat = nan(nSubj,nDelays);
for isubj = 1:nSubj;
    subjnum = subjVec(isubj);
    
    idx = group_data(:,1) == subjnum;
    nTrials = sum(idx);
    
    [b,bint,~,~,stats] = regress(group_data(idx,4), [ones(nTrials,1) group_data(idx,3)]);
    bVec(isubj) = b(2);
    bintVec(isubj,:) = bint(2,:);
    pVec(isubj) = stats(3);
    
    for idelay = 1:nDelays;
        delay = delayVec(idelay);
        idxx = idx & (delay == group_data(:,3));
        
        errorMat(isubj,idelay) = mean(group_data(idxx,4));
    end
      
end

error_m = mean(errorMat);
error_sem = std(errorMat)/sqrt(nSubj);

figure; hold on;
errorbar(delayVec,error_m, error_sem,'k','LineStyle','none');
plot(delayVec,error_m,'ko')
defaultplot

%% PERMUTATION TEST: DELAY TIME

clear all; 

rng(1) % random seed for reproduceability

% load and preprocess data
% ------------------------
subjVec = 4:14;
priorityVec = [0.6 0.3 0.1];

load('exp2_data_trialinfo.mat')
finalerror = sqrt((10.*cosd(group_data(:,13)) - group_data(:,6)).^2 +...
    (10.*sind(group_data(:,13)) - group_data(:,7)).^2);

% columns correspond to [subject, priority, delay, error, circle size]
group_data = [group_data(:,1:3) finalerror group_data(:,4)];
group_data(any(isnan(group_data), 2), :) = [];
   
delayVec = 1000:500:4000;
nPerms = 1000;
nSubj = length(subjVec);
nDelays = length(delayVec);
nPriorities = length(priorityVec);

% combine dleays that are very close to one another
for idelay = 1:nDelays;
    delay = delayVec(idelay);
    
    idx = (delay-100 < group_data(:,3)) & (delay+100 > group_data(:,3));
    group_data(idx,3) = delay;
end

% put each subject's data in terms of z-scores
groupzscoredata = group_data;
for isubj = 1:nSubj;
    subjnum = subjVec(isubj);
    
    idx = group_data(:,1) == subjnum;
    subjdata = group_data(idx,:);
    
    groupzscoredata(idx,4:5) = [zscore(subjdata(:,4)) zscore(subjdata(:,5))];
end

% permute cicle size data for each bin for each subject
nTrialsVec = [];
permdata = nan(size(groupzscoredata,1),nPerms);
for isubj = 1:nSubj; % for each subject
    subjnum = subjVec(isubj);
    idx = groupzscoredata(:,1) == subjnum;
    
    for ipriority = 1:nPriorities % for each priority
        priority = priorityVec(ipriority);
        idxx = idx & (groupzscoredata(:,2) == priority);

        for idelay = 1:nDelays; % for each angle
            delay = delayVec(idelay);
            idxxx = idxx & (groupzscoredata(:,3) == delay); % indices of all trials that fit subject, priority, and angle criteria
            
            nTrials = sum(idxxx); % number of trials that fit criteria
            nTrialsVec = [nTrialsVec nTrials];
            % make a permutation matrix
            permMat = nan(nTrials,nPerms); % permuting these trials nPerms times
            for iperm = 1:nPerms; % for each permutation
                blah = randperm(nTrials);
                while any(blah == 1:nTrials) % implement derangement
                    blah = randperm(nTrials);
                end
                permMat(:,iperm) = blah';
            end
            
            circlesizedata = group_data(idxxx,5); % all wager data for the given subject, priority, and angle
            permdata(idxxx,:) = circlesizedata(permMat); % permuted data 1000 times
        end
    end
    
end

% calculate permuted null distributions
% -------------------------------------

[corrVec, pVec] = deal(nan(1,4*nSubj)); % real correlations and significancee values
[corrMat, pMat] = deal(nan(nPerms,4*nSubj)); % permuted values
for isubj = 1:nSubj; % for each subject
    subjnum = subjVec(isubj);
    idx = groupzscoredata(:,1) == subjnum;

    for ipriority = 1:nPriorities
        priority = priorityVec(ipriority);
         idxx = idx & (groupzscoredata(:,2) == priority);
         
         [r, p] = corrcoef([groupzscoredata(idxx,4:5), permdata(idxx,:)]);
         corrVec(4*(isubj-1) + ipriority) = r(2,1);
         pVec(4*(isubj-1) +ipriority) = p(2,1);
         corrMat(:,4*(isubj-1) +ipriority) = r(3:end,1);
         pMat(:,4*(isubj-1) +ipriority) = p(3:end,1);
    end
end

% print M \pm SEM real and permuated correlations
% -----------------------------------------------

% actual values
corrVec(4:4:44) = [];
data_M = mean(corrVec)
data_SEM = std(corrVec)/sqrt(length(corrVec))

% median of correlation distribution ofpermutations
meds = quantile(corrMat,.5);
meds(4:4:44) = [];
perm_M = mean(meds)
perm_SEM = std(meds)/sqrt(length(meds))


% wilcoxon signed-rank test
% -------------------------
[p,h,stats] = signrank(meds,corrVec)