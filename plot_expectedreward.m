% expected reward plot

%% ======= ONE DIMENSIONAL CASE ===========

% reward function
maxReward = 120;
knobConstant = 33.5738;
slope = 0.4;
rewardFn = @(r)maxReward*exp(-slope*r/knobConstant);

% % plot reward as a function of r
% r = linspace(0,650,50);
% plot(r,rewardFn(r))

% p(S|X)
cushion = 15; % degrees next to cardinal axes in which targets are not presented
p_S = 1/(90-2*cushion); % uniform prior of degree presentation
% S = rand*(90-2*cushion) + cushion; % sample S
% X = S + randn*sigma; % sample X
% p_SgivenX = @(r,sigma)p_S*(normcdf(r,0,sigma) - normcdf(-r,0,sigma)); %p_S*normpdf(X,S,sigma); % posterior
p_SgivenX = @(r,sigma)normcdf(r,0,sigma) - normcdf(-r,0,sigma);

% max expected reward for this sample
nSigma = 50;
sigmaVec = linspace(0,600,nSigma);
for isigma = 1:nSigma;
    sigma = sigmaVec(isigma);
    func = @(r) -rewardFn(r).*p_SgivenX(r,sigma);
    x0 = rand;
    [r(isigma), nll(isigma)] = fminsearch(func,x0);
end

% plot optimal radius as a function of noise
plot(sigmaVec,r,'k-')
defaultplot
ylabel('optimal radius (pixels)')
xlabel('internal memory noise (pixels)')
% title('ONE DIMENSIONAL CASE!')

% figure; 
% plot(sigmaVec,-nll);
% defaultplot;
% ylabel('expected reward')
% xlabel('internal memory noise')

% plot expected utility for given noise as a function of disk size
figure;
sigma = 100;
rVec = 1:500;
EU = rewardFn(rVec).*p_SgivenX(rVec,sigma);
plot(rVec,EU,'k-')
defaultplot;
xlabel('disk size')
ylabel('expected utility')


%% =========== TWO DIMENSIONAL CASE ============
