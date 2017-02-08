function reward = rewardFn(r)
% REWARDFN: reward function of experiment

maxReward = 120; % maxmimum number of points
slope = 0.4;
reward = maxReward*exp(-slope*r);