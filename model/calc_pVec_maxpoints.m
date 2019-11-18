function pVec = calc_pVec_maxpoints(Theta,exppriorityVec)
%CALC_PVEC_MAXPOINTS calculates the proportion allocation that maximizes
%points
% 
%   PVEC = CALC_PVEC_MAXPOINTS(THETA, EXPPRIORITYVEC) returns the optimal 
%     proportion allocation for THETA and EXPPRIORITYVEC. 
% 
%   ========== INPUT VARIABLES ==========
%   
%   THETA: parameter vector [Jbar_total,tau,alpha,beta], where
%       Jbar_total: mean total amount of resources available
%       tau: scale parameter of memory precision gamma distribution
%       alpha: risk preferences for post-decision wager
%       beta: inverse noise temperature on post-decision wager
%
%   EXPPRIORITYVEC: row vector of experimental priority.
%       sum(exppriorityVec) = 1

% ---------------------
%      Aspen H. Yoo
%   aspen.yoo@nyu.edu
% ---------------------

% data stuff
nPriorities = length(exppriorityVec);

% calculate the optimal proportions given the parameters
% calc_ntotalEU = @(x) -(0.6*calc_E_EU([Jbar_total*x(1),tau,alpha,beta]) ...
%     + 0.3*calc_E_EU([Jbar_total*x(2),tau,alpha,beta])...
%     + 0.1*calc_E_EU([Jbar_total*x(3),tau,alpha,beta]));
objfunc = @(x) calc_ntotalEU(Theta,x,exppriorityVec);

% optimization-related variables
Aeq = ones(1,nPriorities);
beq = 1;
[A,b,nonlcon] = deal([]);
options = optimset('Display','none');
lb = 1e-5.*ones(1,nPriorities); 
ub = ones(1,nPriorities);
nStartVals = 15; % tried with different parameters and lowest value showed up 3,5,7,8,10,10 of 10. 

% get starting values
plb = 0.1.*ones(1,nPriorities);
pub = 0.9.*ones(1,nPriorities);
x0 = lhs(nStartVals,nPriorities,plb,pub,[],1e3);
x0 = bsxfun(@rdivide,x0,sum(x0,2)); % normalize

% running optimization nStartVals times
pVec = nan(nStartVals,nPriorities);
nEU = nan(1,nStartVals);
for istartval = 1:nStartVals
    [pVec(istartval,:), nEU(istartval)] = fmincon(objfunc,x0(istartval,:),A,b,Aeq,beq,lb,ub,nonlcon,options);
end

% pick the pVec that corresponded to highest EU
pVec = pVec(nEU == min(nEU),:);
pVec = pVec(1,:); % necessary if multiple entries have nEU == min(nEU)

end

function EU = calc_ntotalEU(Theta,x,exppriorityVec)
Jbar_total = Theta(1);
tau = Theta(2);
alpha = Theta(3);
beta = Theta(4);

EU = 0;
for ipriority = 1:length(exppriorityVec);
    EU = EU - exppriorityVec(ipriority)*calc_E_EU([Jbar_total*x(ipriority),tau,alpha,beta]);
end

end