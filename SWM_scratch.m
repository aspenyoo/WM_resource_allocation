%% understanding Jbar, tau, and sigma

sigma = 5;
Jbar = .04;
tau = .0001;

nSamps = 1e5;
samps = sqrt(1./gamrnd(Jbar/tau,tau,1,nSamps));

hist(1./samps.^2)
% hist(samps)

figure; hist(samps)
mean(samps)
