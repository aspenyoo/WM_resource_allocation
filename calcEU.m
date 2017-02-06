% notes to self
% gneral uncertainty ranges from 1 to 4 ish (up to 5) in visual deg? (disc size?)
% general error 1 to 4 ish (up to 12) in vis deg? (saccade error?)

clear all

% locations
a1 = 10:10:80;
a2 = 100:10:170;
a3 = 190:10:260;
a4 = 280:10:350;
Sdistribution = [a1 a2 a3 a4];
nS = length(Sdistribution);

% parameters
JbarVec = [5 2 1];%1./([2 6 10].^2);       % mean parameter of gamma distribution
tau = 1;%.001;        % scale parameter of gamma distribution
nCond = length(JbarVec); % number of conditions
beta = 1;

% reward function
maxReward = 120;
slope = 0.4;
rewardFn = @(r)maxReward*exp(-slope*r);

nTrials = 1000;
nsamp = 100;
SVec = Sdistribution(randi(nS,nCond,nTrials)); % matrix of stimuli

D = nan(nCond,nTrials); Shat = nan(nCond,nTrials);
rho = nan(1,nCond); pval = nan(1,nCond);
for icond = 1:nCond
    

    Jbar = JbarVec(icond);
    
    % ===== get p(J|Jbar,tau) ======
    % make range small enough
    ind1 = []; ind2 = [];
    blah = 10;
    while isempty(ind2) || isempty(ind2) % make range smaller if nothing is above 1e-10
        blah = blah * 0.9;
        xx = linspace(0,blah,nsamp);
        yy = gampdf(xx,Jbar/tau,tau);
        ind1 = find(yy >= 1e-6,1,'first')-1;
        ind2 = find(yy >= 1e-6,1,'last')+1;
    end
    if ind1 == 0; ind1 = 1; end
    if ind2 == length(yy)+1; ind2 = length(yy); end
    xx = linspace(xx(ind1),xx(ind2),nsamp);
    yy = gampdf(xx,Jbar/tau,tau);
    p_sigma = yy./sum(yy);
    sigval = sqrt(1./xx);
    
    % samples
    sigmaVec = sqrt(1./gamrnd(Jbar/tau,tau,1,nTrials));
    X = SVec(icond,:) + randn(1,nTrials).*sigmaVec;
    Shat(icond,:) = X;
    
    for itrial = 1:nTrials
        sigma = sigmaVec(itrial);
        
        % calculate expected utility across disc size range
        func = @(r) rewardFn(r).*calc_pHit(r,sigma);
        x = linspace(0,10,1e2);
        EUpdf = exp(beta.*func(x)); % pdf of expected utility
        EUpdf = exp(log(EUpdf) - log(sum(EUpdf))); % normalize   
        E_EU = (x.*EUpdf); % expected value of expected utility
        
        
        EUcdf = cumsum(EUpdf);% make cdf
        
        % sample from icdf
        samp = rand;
        idx = find(abs(EUcdf-samp) == min(abs(EUcdf-samp)));
        D(icond,itrial) = x(idx);

    end
    
end