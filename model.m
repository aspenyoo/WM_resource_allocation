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
sigmaVec = [5 7 9];
JbarVec = [5 2 1];%1./(sigmaVec.^2);%1./([2 6 10].^2);       % mean parameter of gamma distribution
tau = 1;%0.008;        % scale parameter of gamma distribution
beta = 1;%15;

nCond = length(JbarVec); % number of conditions

% reward function
maxReward = 120;
slope = 0.4;
rewardFn = @(r)maxReward*exp(-slope*r);

nTrials = 100;
nsamp = 100;
SVec = Sdistribution(ceil(rand(nCond,nTrials).*nS));

colors = ['r';'b';'k'];%aspencolors(nCond,'pastel');%[127 0 0; 247 69 0; 247 148 30]./255; %
figure;
error_sacc = nan(nCond,nTrials);
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
        
        % sample from icdf
        samp = rand;
        idx = find(abs(EUcdf-samp) == min(abs(EUcdf-samp)));
        D(icond,itrial) = x(idx);

    end
    
    % PLOT STUFF!!
    
    % calculte regression and correlation
    %     b = regress(D(icond,:)',[error_sacc(icond,:)' ones(nTrials,1)]);
    
    % plot scatterplot and line
    error_sacc(icond,:) = abs(SVec(icond,:) - Shat(icond,:));
    plot(error_sacc(icond,:),D(icond,:),'o','Color',colors(icond,:));hold on;
    %     blah = [0 10];
    %     plot(blah,b(1).*blah + b(2),colors(icond,:));
    
    [rho(icond), pval(icond)] = corr(error_sacc(icond,:)',D(icond,:)');
end
defaultplot
xlabel('saccadic error');
ylabel('disk size')

rho
pval

% ===== MARGINAL DISTRIBUTIONS ======
nBins = 35;
DMax = max(D(:));
Drange = linspace(0,DMax,nBins);
errorMax = 11;%max(error_sacc(:));
errorrange = linspace(0,errorMax,nBins);

% disk size
figure; hold on;
for icond = 1:nCond;
    [cnts,cntrs] = hist(D(icond,:),Drange);
    plot(cntrs,cnts,'-','Color',colors(icond,:));
end
defaultplot
title('disk size distributions')

% saccadic error
figure; hold on;
for icond = 1:nCond;
    [cnts,cntrs] = hist(error_sacc(icond,:),errorrange);
    plot(cntrs,cnts,'-','Color',colors(icond,:));
end
defaultplot
title('saccadic error distributions')

% ===== BINNED DISK SIZE AND ERROR VAR PLOT =====
nQuants = 4;

figure;
sd_errorsacc = nan(nCond,nQuants);
med_disksize = nan(nCond,nQuants);
for icond = 1:nCond;
   
   % pair disk size and saccade error data.
   dataMat = [D(icond,:); error_sacc(icond,:)]';
   dataMat = sortrows(dataMat);
   
   quantileEnds = linspace(0,nTrials,nQuants+1);
   for iquant = 1:nQuants; % for each bin (by quantile)...
       
       % get appropriate quantile from dataMat
       quantMat = dataMat(quantileEnds(iquant)+1:quantileEnds(iquant+1),:); 
       
       % calculate SD of saccade error for this quantile
       sd_errorsacc(icond,iquant) = std(quantMat(:,2)).^2;
       
       % median of disksizes
       med_disksize(icond,iquant) = median(quantMat(:,1));
       
   end
   
   plot(med_disksize(icond,:),sd_errorsacc(icond,:),'o-','Color',colors(icond,:));
   hold on;
end

defaultplot;
xlabel('disk size')
ylabel('SD of saccade errors')
legend('0.6','0.3','0.1')

% ===== MAIN EFFECT PLOTS (var) =====
sd_disksize = nan(1,nCond);
sd_error = nan(1,nCond);
for icond = 1:nCond;
    sd_disksize(icond) = std(D(icond,:));
    sd_error(icond) = std(error_sacc(icond,:));
end

% condNumVec = [0.6 0.3 0.1];

% saccade error
figure;
plot(sd_error.^2,'k.','MarkerSize',14);
defaultplot;
% xlim([0 0.7])
% ax = gca;
% ax.XTick = condNumVec;
xlabel('priority')
ylabel('SD of saccade errors')

% disksize
figure;
plot(sd_disksize.^2,'k.','MarkerSize',14);
defaultplot
% xlim([0 0.7])
% ax = gca;
% ax.XTick = condNumVec;
xlabel('priority')
ylabel('SD of disk sizes')


