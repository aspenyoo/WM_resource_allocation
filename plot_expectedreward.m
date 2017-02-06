% expected reward plot

%% ======= ONE DIMENSIONAL CASE ===========

% reward function
maxReward = 120;
slope = 0.3;
dva2deg = 1;%pi/18;
rewardFn = @(r)maxReward*exp(-slope*r*dva2deg);

% p(S|X)
cushion = 15; % degrees next to cardinal axes in which targets are not presented
p_S = 1/(90-2*cushion); % uniform prior of degree presentation
p_SgivenX = @(r,sigma)normcdf(r,0,sigma) - normcdf(-r,0,sigma);

% max expected reward for this sample
nSigma = 50;
sigmaVec = linspace(0,10,nSigma);
for isigma = 1:nSigma;
    sigma = sigmaVec(isigma);
    func = @(r) -rewardFn(r).*p_SgivenX(r,sigma);
    x0 = rand;
    [r(isigma), nll(isigma)] = fminsearch(func,x0);
end

% plot optimal radius as a function of noise
plot(sigmaVec,r,'k-')
defaultplot
ylabel('Optimal radius, R* (deg)')
xlabel('Memory uncertainty')
% ax = gca;
%     ax.XTick = [0 5 10];
%     ax.YTick = 0:3;
%     ylim([0 2.5])
    
% plot expected utility for given noise as a function of disk size
figure;
sigma = 1;
rVec = linspace(0,10,100);
EU = rewardFn(rVec).*p_SgivenX(rVec,sigma);
plot(rVec,EU,'k-')
defaultplot;
xlabel('Circular radius (deg)')
ylabel('Expected utility')


%% =========== TWO DIMENSIONAL CASE ============
