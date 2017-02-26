%% understanding Jbar, tau, and sigma

Jbar = 3.5452;
tau = 2;

nSamps = 1e5;
samps = sqrt(1./gamrnd(Jbar/tau,tau,1,nSamps));

% hist(1./samps.^2)
% hist(samps)

figure; hist(samps)
mean(samps)

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


%% plot disc size as a function of euclidean error

isubj = 4;

% nTrials = nan(1,nPriorities);
% for ipriority = 1:nPriorities
%     priority = priorityVec(ipriority);
%     idx = (data_subj == subjnum) & (data_priority == priority);
%     nTrials(ipriority) = sum(idx);
%     
%     data_r = data_radius(idx);
%     data_error = error_distance(idx);
%     
%     subplot(2,2,ipriority)
%     plot(data_error,data_r,'ko')
%     defaultplot
%     
% end

load('cleandata.mat','data')
nPriorities = 3;
model = 2;

Theta = [2.4469    0.5122    0.6620    0.5406    0.2483]; % [2 .2 .2];% load best fit parameters
nTrials = nan(1,3);
for ipriority = 1:nPriorities
    nTrials(ipriority) = size(data{isubj}{ipriority},1);
end

simdata = simulate_data(model,Theta,nTrials);

% plot simulated data on top of real data
for ipriority = 1:nPriorities
    subplot(2,2,ipriority)
    plot(data{isubj}{ipriority}(:,1),data{isubj}{ipriority}(:,2),'ko');
    hold on;
    plot(simdata{ipriority}(:,1),simdata{ipriority}(:,2),'r*');
end

%% plot histogram of error and disc size for each priority

figure;
for ipriority = 1:nPriorities
    
    % plot error
    subplot(3,2,ipriority)
    hist([data{isubj}{ipriority}(:,1) simdata{ipriority}(:,1)])
    
    % plot discsize
    subplot(3,2,3+ipriority)
    hist([data{isubj}{ipriority}(:,2) simdata{ipriority}(:,2)])
end

%% plot median and IQR for error and disc size for each priority

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
