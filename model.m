% notes to self
% gneral uncertainty ranges from 1 to 4 ish (up to 5) in visual deg? (disc size?)
% general error 1 to 4 ish (up to 12) in vis deg? (saccade error?)

close all
clear all
clc

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
SVec = Sdistribution(ceil(rand(nCond,nTrials).*nS));

colors = ['r';'b';'k'];%[127 0 0; 247 69 0; 247 148 30]./255; %
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
%     p_SDgivenXsigma = @(r,sigma) 1./(1+exp(-beta.*rewardFn(r).*p_bringInsidedisc(r,sigma)));
%         p_SgivenX = @(r) sum(p_SDgivenXsigma(r,sigval).*p_sigma);
    for itrial = 1:nTrials;
        sigma = sigmaVec(itrial);
        
        
        % calculate expected utility across disc size range
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
    error_sacc(icond,:) = abs(SVec(icond,:) - Shat(icond,:));
    plot(D(icond,:),error_sacc(icond,:),'o','Color',colors(icond,:));hold on;
    %     blah = [0 10];
    %     plot(blah,b(1).*blah + b(2),colors(icond,:))
    
    % calculate correlation
    [rho(icond), pval(icond)] = corr(error_sacc(icond,:)',D(icond,:)');
    
    %     pause;
    
end
defaultplot
ylabel('saccadic error');
xlabel('disc size')


 
rho
pval

% plot data 
figure;
plot(D(1,:),error_sacc(1,:),'.','Color','k');hold on;


% ===== BINNED disc SIZE AND ERROR SD PLOT =====
nQuants = 4;

figure;
sem_errorsacc = nan(nCond,nQuants); mean_errorsacc = sem_errorsacc;
med_discsize = nan(nCond,nQuants);
for icond = 1:nCond;
   
   % pair disc size and saccade error data.
   dataMat = [D(icond,:); error_sacc(icond,:)]';
   dataMat = sortrows(dataMat);
   
   quantileEnds = linspace(0,nTrials,nQuants+1);
   for iquant = 1:nQuants; % for each bin (by quantile)...
       
       % get appropriate quantile from dataMat
       quantMat = dataMat(quantileEnds(iquant)+1:quantileEnds(iquant+1),:); 
       
       % mean
       mean_errorsacc(icond,iquant) = mean(quantMat(:,2));
       
       % calculate SD of saccade error for this quantile
       sem_errorsacc(icond,iquant) = std(quantMat(:,2))/sqrt(length(quantMat(:,2)));
       
       % median of discsizes
       med_discsize(icond,iquant) = median(quantMat(:,1));
       
   end
   
   errorbar(med_discsize(icond,:),mean_errorsacc(icond,:),sem_errorsacc(icond,:),'Color',colors(icond,:));
   hold on;
end

defaultplot;
xlabel('disc size')
ylabel('mean of saccade errors')
legend('0.6','0.3','0.1')


% ===== MAIN EFFECT PLOTS =====
sem_discsize = nan(1,nCond); mean_discsize = sem_discsize;
sem_error = nan(1,nCond); mean_error = sem_error;
for icond = 1:nCond;
    sem_discsize(icond) = std(D(icond,:))./sqrt(length(D(icond,:)));
    sem_error(icond) = std(error_sacc(icond,:))./sqrt(length(error_sacc(icond,:)));
    mean_discsize(icond) = mean(D(icond,:));
    mean_error(icond) = mean(error_sacc(icond,:));
end

% saccade error
figure;
errorbar(mean_error,sem_error,'k.','MarkerSize',14);
defaultplot
xlabel('priority')
ylabel('mean of saccade errors')

% discsize
figure;
errorbar(mean_discsize,sem_discsize,'k.','MarkerSize',14);
defaultplot
xlabel('priority')
ylabel('mean of disc size')