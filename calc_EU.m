function calc_EU(Jbar,tau,beta)
% calculate expected utility for a given Jbar and tau

% input parameters
Jbar = 4; 
tau = 1;
beta = 1;

% for gamma distribution
nsamp = 100;
xx = linspace(1e-5,10,nsamp);
yy = gampdf(xx,Jbar/tau,tau);

% reward function
maxReward = 120;
slope = 0.4;
rewardFn = @(r)maxReward*exp(-slope*r);

% radius stuff
nsamps_r = 100;
r = linspace(0,3,nsamps_r);

Exp_R = nan(1,nsamp);
for ixx = 1:nsamp;
    J = xx(ixx);
    sigma = sqrt(1/J);

    % function for expected value R
    calc_EU = @(r)calc_p_Hit(r,sigma).*rewardFn(r);
    
    % probability of choosing R.
    p_chooseR = exp(beta.*calc_EU(r));
    p_chooseR = p_chooseR./sum(p_chooseR);
    
    % calculate expected R
    Exp_R(ixx) = calc_EU(r)*p_chooseR(:);
end



% Jbar_total = 10;
% qVec = [0.1 0.3 0.6]; % from low to high
% assert(sum(qVec) == 1)
% 
% JbarVec = qVec*Jbar_total;

