function fit_parameters(testmodel,subjnum,runlist,runmax,truemodel,expnumber,fixparams)
if nargin < 4; runmax = 50; end
if nargin < 5; truemodel = testmodel; end % for model recovery
if isempty(truemodel); truemodel = testmodel; end
if nargin < 6; expnumber = 2; end
if nargin < 7; fixparams = []; end


%
%
% ================= INPUT VARIABLES ==================
% MODEL: 1 (optimal priority placement) or 2 (not optimal)
% SUBJNUM: subject number. 1 - 11
% NSTARTVALS: optimization starting vals
% TESTMODEL: 1 (optimal priority placement) or 2 (not optimal). this is
% used only for model recovery.
% EXPNUMBER: 1 (experiment with just priority manipulation) or 2
% (experiment with disc size response also). 
% 
% -----------------------
%      Aspen H. Yoo
%   aspen.yoo@nyu.edu
%     April 10, 2017


% filepath = ['fits/exp' num2str(expnumber) '_fixedrisk/'];
if isempty(fixparams)
filepath = ['/home/ay963/spatialWM/fits/exp' num2str(expnumber) '/'];
else
filepath = ['/home/ay963/spatialWM/fits/exp' num2str(expnumber) '_fixedrisk/'];
end

if (expnumber == 1) % if nodiscsize experiment (first experiment)
    suffix = '_nodisc';
else
    suffix = [];
end

if (expnumber == 1)
    nSubj = 14;
else
    nSubj = 11;
end

if (subjnum <= nSubj)
    load(['cleandata' suffix '.mat'])
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
pub = [20 5];
logflag = [1 1];
if expnumber == 2 % alpha beta
    lb = [lb 1e-5 1e-5];
    ub = [ub 5 5];
    plb = [plb 0.7 0.5];
    pub = [pub 1.3 1.5];
    logflag = [logflag 0 0];
end
if testmodel == 2 % p_high p_med
    lb = [lb 1e-10 1e-10];
    ub = [ub 1 1];
    plb = [plb 0.3 1e-10];
    pub = [pub 0.7 0.3];
    logflag = [logflag 0 0];
    nonbcon = @(x) sum(x(:,end-1:end),2) > 1;
else
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
    
    try load(filename); catch; ML_parameters = []; nLLVec = []; end
    try blah = runlist_completed*2; catch; runlist_completed = []; end % seeing if runlist_completed exists yet
    
    x0 = x0_list(runlist(irun),:);
%     x0 = plb + rand(1,nParams).*(pub - plb);
%     load(['fits/exp' num2str(expnumber) '/fits_model' num2str(testmodel) '.mat'],'ML_parameters')
%     x0 = ML_parameters(subjnum,:);
    
    fun = @(x) calc_nLL(testmodel,x,subjdata,fixparams);
    
    [x,fval] = bads(fun,x0,lb,ub,plb,pub,nonbcon);

    x(logflag) = exp(x(logflag));
    bfp = nan(1,nParams+size(fixparams,2));
    bfp(freeparamsidx) = x;
    bfp(fixparams(1,:)) = fixparams(2,:);

    ML_parameters = [ML_parameters; bfp];
    nLLVec = [nLLVec fval];
    runlist_completed = [runlist_completed runlist(irun)];
    save(filename,'ML_parameters','nLLVec','runlist_completed')
end


