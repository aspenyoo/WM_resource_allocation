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
JbarVec = [5 2 0.8];%1./([2 6 10].^2);       % mean parameter of gamma distribution
tau = 1;%.001;        % scale parameter of gamma distribution
beta = 1;%15;

nCond = length(JbarVec); % number of conditions

% reward function
maxReward = 120;
slope = 0.4;
rewardFn = @(r)maxReward*exp(-slope*r);

nTrials = 100;
nsamp = 100;
SVec = Sdistribution(ceil(rand(nCond,nTrials).*nS));

colors = aspencolors(nCond,'pastel');
figure;
D = nan(nCond,nTrials); Shat = nan(nCond,nTrials);
rho = nan(1,nCond); pval = nan(1,nCond);
for icond = 1:nCond;
    
    Jbar = JbarVec(icond);
    
    % ===== get p(J|Jbar,tau) ======
    % make range small enough
    ind1 = []; ind2 = [];
    blah = 10;
    while isempty(ind2) || isempty(ind2); % make range smaller if nothing is above 1e-10
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
%     figure; plot(xx,yy); 
    
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
        x = linspace(0,10,1e5);
        EUpdf = exp(beta.*func(x));
        
        EUpdf = EUpdf./sum(EUpdf); % normalize   
        EUcdf = cumsum(EUpdf);% make cdf
        
%         figure;
%         subplot(3,1,1)
%         plot(x,func(x))
%         subplot(3,1,2)
%         plot(x,EUpdf)
%         subplot(3,1,3)
%         plot(x,EUcdf)
        
        % sample from icdf
        samp = rand;
        idx = find(abs(EUcdf-samp) == min(abs(EUcdf-samp)));
        D(icond,itrial) = x(idx);

    end
    
%     figure
    error_sacc = abs(SVec(icond,:) - Shat(icond,:));
    plot(error_sacc,D(icond,:),'o','Color',colors(icond,:));hold on;
    
    % calculate correlation
    [rho(icond), pval(icond)] = corr(error_sacc',D(icond,:)');
    
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


