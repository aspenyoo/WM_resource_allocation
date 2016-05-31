%% understanding Jbar, tau, and sigma

sigma = 9;
Jbar = 1./(sigma.^2); %.04;
tau = .01;

nSamps = 1e5;
samps = sqrt(1./gamrnd(Jbar/tau,tau,1,nSamps));

hist(1./samps.^2)
% hist(samps)

figure; hist(samps)
mean(samps)
