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
filepath = ['fits/exp' num2str(expnumber) '/'];
%  filepath = ['/home/ay963/spatialWM/fits/exp' num2str(expnumber) '/'];
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

[logflag, lb, ub, plb, pub] = loadconstraints(testmodel,expnumber);

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
nParams = length(logflag);
% if ~isempty(nonbcon)
    x0_list = [];
    constantt = 0;
    while size(x0_list,1) < runmax
        rng(0);
        x0_list = lhs(runmax + constantt,nParams,plb,pub,[],1e3);
        countt = nonbcon(testmodel,x0_list);
        x0_list(logical(countt),:) =[];
        constantt = constantt + sum(countt);
    end
% else
%     x0_list = lhs(runmax,nParams,plb,pub,[],1e3);
% end

% optimize for starting values in RUNLIST
for irun = 1:length(runlist)
    runlist(irun)
   
    rng(runlist(irun));
    
    x0 = x0_list(runlist(irun),:);
    fun = @(x) calc_nLL(testmodel,x,subjdata,fixparams);
    
    x = bads(fun,x0,lb,ub,plb,pub,@(y) nonbcon(testmodel,y));
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

