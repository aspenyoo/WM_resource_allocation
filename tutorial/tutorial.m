clear, close all

%% ====================================================================
%            GENERATING MEMORY ESTIMATES WITH THE VP MODEL
% =====================================================================

Jbar = 3;           % mean precision
tau = 0.5;          % shape parameters (related to spread)

%% gamma distribution
% this is the distribution from which memory precision is drawn

nSamps = 5000;       % number of samples

samples = gamrnd(Jbar/tau,tau,[nSamps,1]);
histogram(samples)
xlim([0 10])
xlabel('precision (J)')
ylabel('frequency')

%% how memory noise is drawn

nTrials = 10;

% plot target location 
figure;
plot(0,0,'r.','MarkerSize',24); hold on;
xlabel('x axis'); ylabel('y axis')
axis equal
axis([-2 2 -2 2])

% iteratively plot gaussian SD (sigma) and memory estimates
for itrial = 1:nTrials
    
    % draw a precision from your gamma function
    J = gamrnd(Jbar/tau,tau);
    
    % convert into sigma for gaussian distribution
    sigma = 1/sqrt(J);
    
    % draw a "memory estimate" from your gaussian
    memory = randn(1,2).*sigma;
    
    % plot
    xx = linspace(-sigma,sigma);
    yy = sqrt(sigma^2 - xx.^2);
    plot([xx; xx]',[yy; -yy]',':k');
    pause
    plot(memory(1),memory(2),'ko')
    pause
end

%% simulate a bunch of trials, look at error distributions
% doing the same thing as above, just more efficient...

nTrials = 500;

J = gamrnd(Jbar/tau,tau,[nTrials,1]);               % drawing precisions for each trial
sigma = 1./sqrt(J);                                 % converting to SD of gaussian
memories = bsxfun(@times, randn(nTrials,2),sigma);  % saccade endpoint
error = sqrt(sum(memories.^2,2));                   % euclidean error

% plotting endpoints relative to true location
figure;
plot(0,0,'r.','MarkerSize',24); hold on;
plot(memories(:,1),memories(:,2),'k.')
axis equal
axis([-3 3 -3 3])

figure; 
histogram(error)
xlabel('euclidean error')
ylabel('frequency')

%% ====================================================================
%       NOW LOOKING AT ERROR CLUSTER POINTS FOR MULTIPLE ITEMS 
% =====================================================================
% for now, assuming equal allocation across items

% model parameters
Jbar_total = 10;
tau = 0.5;

% make ups some "experiment information"
nItems = 4;
item_rad = 5;         
item_locs = [1 1; -1 1; -1 -1; 1 -1];
item_locs = item_locs.*item_rad;
nTrials = 100;

%% simulate data for all items

% precision for each item
Jbar = Jbar_total./nItems; % this is where the equal allocation comes in

memories = cell(1,nItems);
for iitem = 1:nItems;
    item_loc = item_locs(iitem,:);
    
    % generate memories
    % 00000000000000 !! YOU GOTTA MAKE THIS FUNCTION !! 00000000000
    memories{iitem} = simulate_memories(Jbar,tau,item_loc,nTrials);
    % 00000000000000000000000000000000000000000000000000000000000000
end

%% plot the simulated data

figure; hold on
axis equal
axis(10*[-1 1 -1 1])

% plot memories
for iitem = 1:nItems
    plot(memories{iitem}(:,1),memories{iitem}(:,2),'k.')
end

% plot true item locations
plot(item_locs(:,1),item_locs(:,2),'r.','MarkerSize',24);

%% ====================================================================
%                    UNEQUAL ALLOCATION MODEL 
% =====================================================================
colorMat = [1 0 0; 0 0 1; 0 0 0];

%% looking at gamma distributions based on allocation

allocationVec = [0.6 0.3 0.1];
xx = linspace(0,10,100);

figure; hold on
for ipriority = 1:length(allocationVec)
    Jbar = Jbar_total*allocationVec(ipriority);
    
    plot(xx,gampdf(xx,Jbar/tau,tau),'Color',colorMat(ipriority,:))
end
xlabel('precision (J)')
ylabel('proportion')

%% simulate memories based on some allocation vector

% total number of trials
nTrialsTotal = 1000;

% actual priority values (include only nonzero priorities)
expPriorityVec = [0.6 0.3 0.1];

% proportion of resource allocated to items
allocationVec = []; % <--- 0000 FILL THIS IN 00000

% 000000000000000000 !! FILL THE STUFF IN BELOW !! 00000000000000
% simulate memory data
nPriorities = length(expPriorityVec);    % number of nonzero priority items
memories = cell(1,nPriorities);
for ipriority = 1:nPriorities    % for each item...
    item_loc = [];      % item location
    priority = [];      % item priority
    p = [];             % proportion allocated to item
    Jbar = [];          % item precision
    N = [];             % number of trials for current item
    
    % generate memory
    memories{ipriority} = simulate_memories(Jbar,tau,item_loc,N);
end
% 00000000000000000000000000000000000000000000000000000000000000

%% plot the simulated data

figure; hold on
axis equal
axis(10*[-1 1 -1 1])

% plot memories
for ipriority = 1:nPriorities
    plot(memories{ipriority}(:,1),memories{ipriority}(:,2),'k.')
end

% plot true item locations
plot(item_locs(:,1),item_locs(:,2),'r.','MarkerSize',24);

%% ====================================================================
%              ESTIMATING PARAMETERS: PROPORTIONAL MODEL 
% =====================================================================
% the proportional model has two parameters: Jbar_total and tau

%% simplifying error

% 0000000000000000000000 !! FILL IN !! 0000000000000000000000
% should be a cell of length nPriorities
error = simulate_data(model,expnumber,Theta,nTrials,expPriorityVec);
% 00000000000000000000000000000000000000000000000000000000000000

% plot
xlims = linspace(0,10,16); % x values for histogram
figure; hold on;
for ipriority = 1:nPriorities
    datacounts = hist(error{ipriority},xlims);
    plot(xlims,datacounts./sum(datacounts),'Color',colorMat(ipriority,:));
end

%% calculating log likelihood of a parameter combination

% load data
load('exp1_cleandata.mat')
subjnum = 5;                    % subject number
data = data{subjnum};   

model = 'proportional';
Theta = [Jbar_total tau];
fixparams = [];

% calculate -LL
nLL = calc_nLL(model,Theta,data,expPriorityVec,fixparams)

%% fitting parameters       

% model fitting stuff
model = 'proportional';             % model name
exppriorityVec = [0.6 0.3 0.1];            % experimental priority vector
runlist = 1;                    % ignore. which idxs of total runs for current model/data
runmax = 20;                    % ignore. number of runs per model/data
fixparams = [];                 % fixed parameters, ignore for now

% fit parameter
[ML_parameters, nLLVec] = fit_parameters(model,data,exppriorityVec,runlist,runmax,fixparams);

%% plotting model fits

% plot data
figure; hold on;
for ipriority = 1:nPriorities
    datacounts = hist(data{ipriority},xlims);
    plot(xlims,datacounts./sum(datacounts),'Color',colorMat(ipriority,:));
end

% get model prediction
expnumber = 1;
error = simulate_data(model,expnumber,ML_parameters,nTrials,expPriorityVec);
error = error{1};

% plot model prediction
for ipriority = 1:nPriorities
    datacounts = hist(error{ipriority},xlims);
    plot(xlims,datacounts./sum(datacounts),':','Color',colorMat(ipriority,:));
end

%% ====================================================================
%                    CALCULATING EXPECTED LOSS
% =====================================================================

clear all

%% for a single item

% getting parameters
Jbar = 3;
tau = 1;
gamma = 1;

% E[d^gamma] = \int d^blah \int p(d|J)p(J) dd dJ
% E[d^gamma] = \int d^blah \int rayleigh(d|1/J) gamma(J|Jbar,tau) dd dJ

% get x-axis of J and distance d
[JVec,dVec] = loadvar('JVec',{Jbar,tau},'rVec');
JVec = JVec';

% get Jpdf: p(J|Jbar,tau). 500 x 1 vector
Jpdf = gampdf(JVec,Jbar/tau,tau); % probability of that J value
Jpdf = Jpdf./sum(Jpdf);

% get dpdf given J: p(d|1/J). will be nJs (500) x nds (500) vector
d_given_J = bsxfun(@(x,y) x.*y.*exp(-x^2.*y/2),dVec,JVec);

% get dpdf (marginalize over J)
dpdf = bsxfunandsum(@times,d_given_J,Jpdf);
dpdf = dpdf./sum(dpdf);

% \int d^gamma p(d) dd
expectederror = sum((dVec.^gamma) .* dpdf)

%% across items

% getting parameters
Jbar_total = 10;
tau = 1;
gamma = 1;

priorityVec = [0.6 0.3 0.1];
allocatedpriorityVec = priorityVec;
nPriorities = length(allocatedpriorityVec);

% E[d^gamma] = \int d^blah \int p(d|J)p(J) dd dJ
% E[d^gamma] = \int d^blah \int rayleigh(d|1/J) gamma(J|Jbar,tau) dd dJ

expectederror = 0;
for ipriority = 1:nPriorities
    Jbar = Jbar_total*allocatedpriorityVec(ipriority);
    
    % get x-axis of J and distance d 
    [JVec,dVec] = loadvar('JVec',{Jbar,tau},'rVec');
    JVec = JVec';
    
    % get Jpdf: p(J|Jbar,tau). 500 x 1 vector
    Jpdf = gampdf(JVec,Jbar/tau,tau); % probability of that J value
    Jpdf = Jpdf./sum(Jpdf);
    
    % get dpdf given J: p(d|1/J). will be nJs (500) x nds (500) vector
    d_given_J = bsxfun(@(x,y) x.*y.*exp(-x^2.*y/2),dVec,JVec);
    
    % get dpdf (marginalize over J)
    dpdf = bsxfunandsum(@times,d_given_J,Jpdf);
    dpdf = dpdf./sum(dpdf);
    
    % \int d^blah p(d) dd
    expectederror = expectederror + priorityVec(ipriority).*sum((dVec.^gamma) .* dpdf);
    
end

% see function expectederror = calc_expectederror(Theta,allocatedpriorityVec)
% expectederror = calc_expectederror_analytical(Theta,allocatedpriorityVec,exppriorityVec)

%% ====================================================================
%                    MINIMIZING "ERROR"
% =====================================================================

[pVec, fval] = calc_pVec_minerror(Theta,exppriorityVec)

%% ====================================================================
%                    FITTING MINIMIZING ERROR MODEL
% =====================================================================

%% ====================================================================
%                    CALCULATING EXPECTED UTILITY (EU) 
% =====================================================================

%% EU OF A SINGLE TRIAL


%% EU ACROSS THE EXPERIMENT (EXPECTED EXPECTED UTILITY??)


%% ====================================================================
%                 CALCULATING OPTIMAL CIRCLE SIZE
% =====================================================================


%% ====================================================================
%                       MAXIMIZING POINTS 
% =====================================================================

%% ====================================================================
%                 FITTING MAXIMIZING POINTS MODEL
% =====================================================================
