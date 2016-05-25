% notes to self
% gneral uncertainty ranges from 1 to 4 ish (up to 5) in visual deg? (disk size?)
% general error 1 to 4 ish (up to 12) in vis deg? (saccade error?)


% locations
a1 = 10:10:80;
a2 = 100:10:170;
a3 = 190:10:260;
a4 = 280:10:350;
Sdistribution = [a1 a2 a3 a4];
nS = length(Sdistribution);

% parameters
JbarVec = [25 20 10 5];       % shape parameter of gamma distribution
tau = 3;        % scale parameter of gamma distribution

nCond = length(JbarVec); % number of conditions
% % plot gamma distribution
% samps = gamrnd(Jbar,tau,10000,1);
% hist(samps,50)

% reward function
maxReward = 120;
knobConstant = 33.5738;
slope = 0.4;
rewardFn = @(r)maxReward*exp(-slope*r/knobConstant);

nTrials = 100;
SVec = Sdistribution(ceil(rand(nCond,nTrials).*nS));

colors = [.6 .2 .6; .7 .7 .1; .8 .3 0; 0 .4 .7]; %aspencolors(nCond,'pastel');
figure; 
r = nan(nCond,nTrials); Shat = nan(nCond,nTrials);
rho = nan(1,nCond); pval = nan(1,nCond);
for icond = 1:nCond;
    
    Jbar = JbarVec(icond);
    
    % samples
    sigmaVec = sqrt(1./gamrnd(Jbar/tau,tau,1,nTrials));
    X = SVec(icond,:) + randn(1,nTrials).*sigmaVec;
    Shat(icond,:) = X;
    
    p_SgivenX = @(r,sigma)normcdf(r,0,sigma) - normcdf(-r,0,sigma);
    for itrial = 1:nTrials;
        sigma = sigmaVec(itrial);
        
        % max expected reward for this sample
        func = @(r) -rewardFn(r).*p_SgivenX(r,sigma);
        x0 = rand;
        r(icond,itrial) = fminsearch(func,x0);
    end
    
    error_sacc = abs(SVec(icond,:) - Shat(icond,:));
    plot(error_sacc,r(icond,:),'o','Color',colors(icond,:));hold on;
    
    % calculate correlation
    [rho(icond), pval(icond)] = corr(error_sacc',r(icond,:)');
    
%     pause;
    
end
defaultplot
xlabel('saccadic error');
ylabel('disk size')


rho
pval

% MARGINAL DISTRIBUTIONS

% % disk size ("general uncertainty")
% figure; hold on;
% for icond = 1:nCond;
%     [cnts,cntrs] = hist(r(icond,:),10);
%     plot(cntrs,cnts,'-','Color',colors(icond,:));
% end
% defaultplot
% title('disk size distributions')
% 
% % saccadic error ("general error")
% figure; hold on;
% for icond = 1:nCond;
%     [cnts,cntrs] = hist(abs(SVec(icond,:) - Shat(icond,:)),10);
%     plot(cntrs,cnts,'-','Color',colors(icond,:));
% end
% defaultplot
% title('saccadic error distributions')


