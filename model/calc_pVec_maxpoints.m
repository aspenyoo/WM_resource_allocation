function pVec = calc_pVec_maxpoints(Theta)
%calc_pVec_maxpoints calculates the proportion allocation that maximizes
%points
% 
%   PVEC = CALC_PVEC_MAXPOINTS(THETA) returns the optimal proportion
%   allocation for THETA = [Jbar_total, tau, alpha, beta], where
%       Jbar_total: mean total amount of resources available
%       tau: scale parameter of memory precision gamma distribution
%       alpha: risk preferences for post-decision wager
%       beta: inverse noise temperature on post-decision wager

% ---------------------
%      Aspen H. Yoo
%   aspen.yoo@nyu.edu
% ---------------------


Jbar_total = Theta(1);
tau = Theta(2);
alpha = Theta(3);w
beta = Theta(4);

% data stuff
priorityVec = [0.6 0.3 0.1];
nPriorities = length(priorityVec);

% calculate the optimal proportions given the parameters
calc_ntotalEU = @(x) -(0.6*calc_E_EU([Jbar_total*x(1),tau,alpha,beta]) ...
    + 0.3*calc_E_EU([Jbar_total*x(2),tau,alpha,beta])...
    + 0.1*calc_E_EU([Jbar_total*x(3),tau,alpha,beta]));

% optimization-related variables
Aeq = [1 1 1];
beq = 1;
[A,b,nonlcon] = deal([]);
options = optimset('Display','none');
lb = [1e-5 1e-5 1e-5];
ub = [1 1 1];
nStartVals = 10; % tried with different parameters and lowest value showed up 3,5,7,8,10,10 of 10. 

% running optimization nStartVals times
pVec = nan(nStartVals,3);
nEU = nan(1,nStartVals);
for istartval = 1:nStartVals
    [pVec(istartval,:), nEU(istartval)] = fmincon(calc_ntotalEU,rand(1,3),A,b,Aeq,beq,lb,ub,nonlcon,options);
end

% pick the pVec that corresponded to highest EU
pVec = pVec(nEU == min(nEU),:);
pVec = pVec(1,:); % necessary if multiple entries have nEU == min(nEU)

