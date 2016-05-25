function nll = loglikelihood(data,parameterVec)

% parameters
JbarVec = parameterVec(1:end-2);
tau = parameterVec(end-1);
beta = parameterVec(end);
nCond = length(JbarVec);
if length(data)~=nCond; error('number of parameters does not match length of data cell!'); end
diskSizeVec = 1:120;

% reward function
maxReward = 120;
knobConstant = 33.5738; % note that data collected after 4/25/16 will have a different knob constant!
slope = 0.4;
rewardFn = @(r)maxReward*exp(-slope*r/knobConstant);

% calculate log likelihood
nll_disksize = nan(1,nCond); nll_sacc = nan(1,nCond);
for icond = 1:nCond;
    
    S = data{icond}(:,1)';              % stimulus location (polar angle)
    Shat = data{icond}(:,2)';           % subjects' saccaded to location (polar angle)
    disksize = data{icond}(:,3)';       % subject's reported disk size (radius in pixels)
    nTrials = size(S,1);
    
    Jbar = JbarVec(icond);
    
    % samples
    sigma = gamrnd(Jbar/tau,tau,1,nTrials);
    X = S + randn(1,nTrials).*sigma;
    
    p_SgivenX = @(r,sigma)normcdf(r,0,sigma) - normcdf(-r,0,sigma);
    
    % p(r|X,sigma)
    EU = @(r,sig) rewardFn(r).*p_SgivenX(r,sig); % expected utility
    cellsigma = num2cell(sigma);
    denom = cellfun(@(x) sum(exp(beta.*EU(diskSizeVec,x*ones(1,length(diskSizeVec))))), cellsigma,'UniformOutput',false);
    denom = cell2mat(denom);
    nll_disksize(icond) = -sum(beta.*EU(disksize,sigma) - log(denom));
    
    % p(Shat|X,sigma)
    nll_sacc(icond) = -sum(log(normpdf(Shat,X,sigma)));

    
    nll = sum(nll_disksize) + sum(nll_sacc);

end