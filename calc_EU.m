function E_EU = calc_EU(Jbar,tau,beta)
% CALC_EU calculates expected value of the expected utility for a given 
% Jbar, tau, and beta. 
%
% aspen yoo -- 02.07.2016

% % input parameters
% Jbar = 15; 
% tau = 1;
% beta = 1;

% for gamma distribution
JVec = loadvar('JVec');
nJs = length(JVec);

Jpdf = gampdf(JVec,Jbar/tau,tau);
Jpdf = Jpdf./sum(Jpdf);

% reward function
maxReward = 120;
slope = 0.4;
rewardFn = @(r)maxReward*exp(-slope*r);

% radius stuff
nRs = 100;
rVec = linspace(0,3,nRs); % ASPEN: make sure this range is reasonable

% EU(r,J) for each combination of r and J in rVec and jVec
EU = bsxfun(@times, calc_p_Hit(rVec,JVec), rewardFn(rVec)');

% p(r)
p_chooseR = exp(beta.*EU);
p_chooseR = bsxfun(@rdivide, p_chooseR, sum(p_chooseR)); % normalize across rVec

% EU(J) = \int EU(r,J) p(r) dr
EU_J = sum(EU.*p_chooseR);

% EU(Jbar) = \int EU(J) p(J) dJ
E_EU = Jpdf*EU_J';
