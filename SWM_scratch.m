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

% priority and subject numbers
priorityVec = sort(unique(data_priority),'descend'); % priority condition
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
        
        data{isubj}{ipriority} = nan(nTrials,2);
        data{isubj}{ipriority}(:,1) = error_distance(idx);
        data{isubj}{ipriority}(:,2) = data_radius(idx);
    end
end

save('cleandata.mat','data')

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%    LOOK AT FITTED PARAMETERS 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

clear all

load('cleandata.mat')

for isubj = [1 2 4:11];
isubj

% load fits
filename = ['fits_optimal_subj' num2str(isubj) '.mat'];
load(filename)

bfp = ML_parameters(nLLVec == min(nLLVec),:);
logflag = logical([1 1 0]);
logbfp = bfp;
logbfp(logflag) = log(logbfp(logflag));

pVec = calculate_optimal_pVec(bfp)
min(nLLVec)
nLL = calc_nLL(logbfp,data{isubj})
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % 
%       GET A PLOT OF DATA AND MODEL 
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
clear all
load('cleandata.mat')
nPriorities = 3;

isubj = 1;
subjdata = data{isubj};
for ipriority = 1:nPriorities;
    priority = priorityVec(ipriority);

    
end


