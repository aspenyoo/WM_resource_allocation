%% PSEUDO CODE FOR MODEL

J = 0.06; %
sigma = sqrt(1/J);
beta = 1;

% reward function
maxReward = 120;
slope = 0.4;
rewardFn = @(r)maxReward*exp(-slope*r);

%% doing stuff
nsamps_r = 100;
r = linspace(0,3,nsamps_r);

% function for expected value R
calc_EU = @(r)calc_p_Hit(r,sigma).*rewardFn(r);

% probability of choosing R. 
p_chooseR = exp(beta.*calc_EU(r));
p_chooseR = p_chooseR./sum(p_chooseR);

% calculate expected R
Exp_R = calc_EU(r)*p_chooseR(:);


