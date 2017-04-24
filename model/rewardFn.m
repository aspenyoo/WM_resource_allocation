function reward = rewardFn(r,alpha)
% REWARDFN: reward function of experiment
% 
% ========= INPUT VARIABLES =======
% R: radius of confidence disc size. can be a vector
% ALPHA: risk-averse risk-seeking parameter

maxReward = 120; % maxmimum number of points
slope = 0.4;
reward = maxReward*exp(-slope*r*alpha);