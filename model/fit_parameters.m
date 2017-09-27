function fit_parameters(testmodel,subjnum,runlist,runmax,truemodel,expnumber,fixparams)
if nargin < 4; runmax = 50; end
if nargin < 5; truemodel = testmodel; end % for model recovery
if isempty(truemodel); truemodel = testmodel; end
if nargin < 6; expnumber = 2; end
if nargin < 7; fixparams = []; end


% ================= INPUT VARIABLES ==================
% MODEL: 1 (optimal priority placement) or 2 (not optimal) or 3 (fixed)
% SUBJNUM: subject number. 1 - 11
% NSTARTVALS: optimization starting vals
% TESTMODEL: 1 (optimal priority placement) or 2 (not optimal) or 3 (fixed). 
% this is used only for model recovery.
% EXPNUMBER: 1 (experiment with just priority manipulation) or 2
% (experiment with disc size response also). 
% 
% -----------------------
%      Aspen H. Yoo
%   aspen.yoo@nyu.edu
%     April 10, 2017


% filepath = ['fits/exp' num2str(expnumber) '_fixedrisk/'];
if isempty(fixparams)
% filepath = ['fits/exp' num2str(expnumber) '/'];
 filepath = ['/home/ay963/spatialWM/fits/exp' num2str(expnumber) '/'];
else
filepath = ['/home/ay963/spatialWM/fits/exp' num2str(expnumber) '_fixedrisk/'];
end

if (expnumber == 1)
    nSubj = 14;
else
    nSubj = 11;
end

if (subjnum <= nSubj)
    load(['exp' num2str(expnumber) '_cleandata.mat'])
    subjdata = data{subjnum};
    filename = [filepath 'fits_model' num2str(testmodel) '_subj' num2str(subjnum) '.mat'];
else
    load([filepath 'simdata_model' num2str(truemodel) '.mat'],'simdata')
    subjdata = simdata{subjnum - nSubj};
    if strcmp(truemodel,testmodel) % if same model
    filename = [filepath 'paramrecov_model' num2str(testmodel) '_subj' num2str(subjnum-nSubj) '.mat'];
    else
        filename = [filepath 'modelrecov_truemodel' num2str(truemodel) '_testmodel' num2str(testmodel) '_subj' num2str(subjnum-nSubj) '.mat'];
    end
end

% lower and upper bounds, logflags
lb = [1e-5 1e-3]; % Jbar_total, tau
ub = [50 10];
plb = [0.5 0.01];
pub = [10 1];
logflag = [1 1];
if expnumber == 2 % alpha beta
    lb = [lb 1e-5 1e-5];
    ub = [ub 5 5];
    plb = [plb 0.7 0.5];
    pub = [pub 1.3 1.5];
    logflag = [logflag 0 0];
end
switch testmodel
    case 2 % define p_high p_med for flexible model
        lb = [lb 1e-10 1e-10];
        ub = [ub 1 1];
        plb = [plb 0.3 1e-10];
        pub = [pub 0.7 0.3];
        logflag = [logflag 0 0];
        nonbcon = @(x) sum(x(:,end-1:end),2) >= 1; % violates if p_high + p_med >= 1
    case 4
        lb = [lb 1e-10];
        ub = [ub 10];
        plb = [plb 1e-3];
        pub = [pub 1];
        logflag = [logflag 1];
        nonbcon = @model4nonbcon; % violates if Jbar/tau - psi/2 <=0 
    otherwise
        nonbcon = [];
end
logflag = logical(logflag); 
lb(logflag) = log(lb(logflag)); 
ub(logflag) = log(ub(logflag)); 
plb(logflag) = log(plb(logflag)); 
pub(logflag) = log(pub(logflag)); 

if ~(isempty(fixparams))
    % vector of non-fixed parameter indices
    nParams= length(logflag);
    freeparamsidx = 1:nParams;
    freeparamsidx(fixparams(1,:)) = [];
    
    logflag(fixparams(1,:)) = [];
    lb(fixparams(1,:)) = [];
    ub(fixparams(1,:)) = [];
    plb(fixparams(1,:)) = [];
    pub(fixparams(1,:)) = [];
end

% create list of all x0s
rng(0);
nParams = length(logflag);
x0_list = lhs(runmax,nParams,plb,pub,[],1e3);

% optimize for starting values in RUNLIST
for irun = 1:length(runlist)
    runlist(irun)
   
    rng(runlist(irun));
    
    x0 = x0_list(runlist(irun),:);
    fun = @(x) calc_nLL(testmodel,x,subjdata,fixparams);
    
    x = bads(fun,x0,lb,ub,plb,pub,nonbcon);
    fval = fun(x);
    
    x(logflag) = exp(x(logflag));
    if isempty(fixparams)
        bfp = x;
    else
        bfp = nan(1,nParams+size(fixparams,2));
        bfp(freeparamsidx) = x;
        bfp(fixparams(1,:)) = fixparams(2,:);
    end

    try load(filename); catch; ML_parameters = []; nLLVec = []; end
    try runlist_completed*2; catch; runlist_completed = []; end % seeing if runlist_completed exists yet
    
    ML_parameters = [ML_parameters; bfp];
    nLLVec = [nLLVec fval];
    runlist_completed = [runlist_completed runlist(irun)];
    save(filename,'ML_parameters','nLLVec','runlist_completed')
end
end

function countt = model4nonbcon(x)
    countt = 0;
    countt = countt + (exp(x(:,1))./exp(x(:,2))) <= (exp(x(:,end))./2); % k > psi/2
    countt = countt + (exp(x(:,end)).*3 > exp(x(:,1))); % Jbar_total > psi*3
    countt = countt + exp(x(:,1)) <= 3.*exp(x(:,2)); % Jbar_total > 3*tau
end

