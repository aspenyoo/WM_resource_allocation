% expected reward plot

%% ======= ONE DIMENSIONAL CASE ===========

% reward function
maxReward = 120;
slope = 0.4;
rewardFn = @(r)maxReward*exp(-slope*r);

% p(S|X)
cushion = 15; % degrees next to cardinal axes in which targets are not presented
p_S = 1/(90-2*cushion); % uniform prior of degree presentation
p_SgivenX = @(r,sigma)normcdf(r,0,sigma) - normcdf(-r,0,sigma);

% max expected reward for this sample
nSigma = 50;
sigmaVec = linspace(0,20,nSigma);
for isigma = 1:nSigma;
    sigma = sigmaVec(isigma);
    func = @(r) -rewardFn(r).*p_SgivenX(r,sigma);
    x0 = rand;
    [r(isigma), nll(isigma)] = fminsearch(func,x0);
end

% plot optimal radius as a function of noise
plot(sigmaVec,r,'k-')
defaultplot
ylabel('optimal radius, D* (dva)')
xlabel('uncertainty')
ax = gca;
    ax.XTick = [0 10 20];
    
% plot expected utility for given noise as a function of disk size
figure;
sigma = 6;
rVec = linspace(0,25,100);
EU = rewardFn(rVec).*p_SgivenX(rVec,sigma);
plot(rVec,EU,'k-')
defaultplot;
xlabel('disk size, D (dva)')
ylabel('expected utility')


%% =========== TWO DIMENSIONAL CASE ============
