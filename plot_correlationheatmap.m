% notes to self
% gneral uncertainty ranges from 1 to 4 ish (up to 5) in visual deg? (disk size?)
% general error 1 to 4 ish (up to 12) in vis deg? (saccade error?)


% locations
a1 = 10:10:80;
a2 = 100:10:170;
a3 = 190:10:260;
a4 = 280:10:350;
Sdistribution = [a1 a2 a3 a4]';
nS = length(Sdistribution);

% parameters
nJbar = 20;
nTau = 20;
JbarVec = exp(linspace(-3,3.5,nJbar));       % shape parameter of gamma distribution
tauVec =  exp(linspace(-3,3.5,nTau));         % scale parameter of gamma distribution

% reward function
maxReward = 120;
knobConstant = 33.5738;
slope = 0.4;
rewardFn = @(r)maxReward*exp(-slope*r/knobConstant);

nTrials = 500;

mean_saccerror = nan(nJbar,nTau); mean_disksize  = mean_saccerror;
median_saccerror = nan(nJbar,nTau); median_disksize  = median_saccerror;
mode_saccerror = nan(nJbar,nTau); mode_disksize  = mode_saccerror;
corrcoeff = nan(nJbar,nTau); pval = nan(nJbar,nTau);
for ijbar = 1:nJbar;

    ijbar
    Jbar = JbarVec(ijbar);
    
    for itau = 1:nTau;
        tau = tauVec(itau);
        
        % samples
        SVec = Sdistribution(ceil(rand(nTrials,1).*nS));
        sigmaVec = sqrt(1./gamrnd(Jbar/tau,tau,nTrials,1));
        X = SVec + randn(nTrials,1).*sigmaVec;
        Shat = X;
        
        p_SgivenX = @(r,sigma)normcdf(r,0,sigma) - normcdf(-r,0,sigma);
        r = nan(nTrials,1);
        for itrial = 1:nTrials;
            sigma = sigmaVec(itrial);
            
            % max expected reward for this sample
            func = @(r) -rewardFn(r).*p_SgivenX(r,sigma);
            x0 = rand;
            r(itrial) = fminsearch(func,x0);
        end
        
        error_sacc = abs(SVec - Shat);
        
        % calculate correlation
        [corrcoeff(ijbar,itau), pval(ijbar,itau)] = corr(error_sacc,r);
        
        % measures of central tendencies
        mean_disksize(ijbar,itau) = mean(r);
        mean_saccerror(ijbar,itau) = mean(error_sacc);
        median_disksize(ijbar,itau) = median(r);
        median_saccerror(ijbar,itau) = median(error_sacc);
        
        mode_disksize(ijbar,itau) = mean(mode(ceil(r)));
        mode_saccerror(ijbar,itau) = mean(mode(ceil(error_sacc)));
    end
    
end

imagesc(JbarVec,tauVec,corrcoeff); colormap('parula'); defaultplot;
ylabel('Jbar')
xlabel('tau')
title('Pearson r')

figure;
imagesc(JbarVec,tauVec,pval); colormap('parula'); defaultplot
ylabel('Jbar')
xlabel('tau')
title('pvalue of correlation coefficient')

figure; 
imagesc(JbarVec(4:end),tauVec(1:end-3),mode_saccerror(4:end,1:end-3)); colormap('parula'); defaultplot
ylabel('Jbar')
xlabel('tau')
title('mode saccadic error')

figure; 
imagesc(JbarVec,tauVec,median_disksize); colormap('parula'); defaultplot
ylabel('Jbar')
xlabel('tau')
title('median disk size')