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
beta = 1;

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
nsamp = 100;
SVec = Sdistribution(ceil(rand(nCond,nTrials).*nS));

colors = [.6 .2 .6; .7 .7 .1; .8 .3 0; 0 .4 .7]; %aspencolors(nCond,'pastel');
figure;
r = nan(nCond,nTrials); Shat = nan(nCond,nTrials);
rho = nan(1,nCond); pval = nan(1,nCond);
for icond = 1:nCond;
    
    Jbar = JbarVec(icond);
    
    % get p(J|Jbar,tau)
    xx = linspace(0,200,nsamp);
    yy = gampdf(xx,Jbar/tau,tau);
    ind1 = find(yy >= 1e-10,1,'first') - 1;
    ind2 = find(yy >= 1e-10,1,'last') + 1;
    xx = linspace(xx(ind1),xx(ind2),nsamp);
    yy = gampdf(xx,Jbar/tau,tau);
    p_sigma = yy./sum(yy);
    sigval = sqrt(1./xx);
    
    % samples
    sigmaVec = sqrt(1./gamrnd(Jbar/tau,tau,1,nTrials));
    X = SVec(icond,:) + randn(1,nTrials).*sigmaVec;
    Shat(icond,:) = X;
    
    p_SgivenX = @(r,sigma) (normcdf(r,0,sigma) - normcdf(-r,0,sigma));
%     p_SDgivenXsigma = @(r,sigma) 1./(1+exp(-beta.*rewardFn(r).*p_bringInsideDisk(r,sigma)));
%         p_SgivenX = @(r) sum(p_SDgivenXsigma(r,sigval).*p_sigma);
    for itrial = 1:nTrials;
        sigma = sigmaVec(itrial);
        
        
        % calculate expected utility across disk size range
        func = @(r) rewardFn(r).*p_SgivenX(r,sigma);
        xx = linspace(0,25,1e5);
        EUpdf = exp(beta.*func(xx));
        
        EUpdf = EUpdf./sum(EUpdf); % normalize   
        EUcdf = cumsum(EUpdf);% make cdf
        
        % sample from icdf
        samp = rand;
        idx = find(abs(EUcdf-samp) == min(abs(EUcdf-samp)));
        r(icond,itrial) = xx(idx);

%         % max expected reward for this sample
%         %         
%         func = @(r) -p_SDgivenXsigma(r,sigma);
%         x0 = rand;
%         r(icond,itrial) = fminsearch(func,x0);
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


