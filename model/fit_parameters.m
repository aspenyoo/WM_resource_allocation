function [ML_parameters, nLLVec, runlist_completed] = fit_parameters(model,data,exppriorityVec,runlist,runmax,fixparams)
%FIT_PARAMETERS
% 
%   FIT_PARAMETERS(MODEL, DATA, EXPPRIORITYVEC) estimates and saves
%     maximum-likelihood parameter estimate for MODEL, DATA, EXPPRIORITYVEC, and
%     RUNLIST.
% 
%   FIT_PARAMETERS(MODEL, DATA, EXPPRIORITYVEC, RUNLIST, RUNMAX) estimates and saves
%     maximum-likelihood parameter estimate for TESTMODEL, SUBJNUM,
%     EXPPRIORITYVEC, and RUNLIST.
% 
%   ===== INPUT VARIABLES =====
% 
%   MODEL: 'max_points', 'flexible', 'proportional', or 'min_error'
% 
%   DATA: struct of length nPriorities, with column vector(s) of data (in decreasing order of priority) 
%
%   EXPPRIORITYVEC: vector of priority values (in decreasing order of
%   priority) e.g., [0.6 0.3 0.1]
% 
%   RUNLIST: scalar or vector containing indexes 1 to RUNMAX.
% 
%   RUNMAX: total number of optimizations you plan to start for a given
%     TESTMODEL and SUBJNUM. used in seed, necessary for replicability. 
% 
%   FIXPARAMS: fixed parameters. 2 x (number of fixed parameters), where
%     first row corresponds to the index and second row corresponds to the
%     value of the fixed parameter. 

% -----------------------
%      Aspen H. Yoo
%   aspen.yoo@nyu.edu
% -----------------------

if nargin < 4; runlist = 1; end
if nargin < 5; runmax = 50; end
if nargin < 6; fixparams = []; end

% exppriorityVec = [0.6 0.3 0.1];
expnumber = size(data{1},2);
assert(length(data) == length(exppriorityVec));
[logflag, lb, ub, plb, pub] = loadconstraints(model,expnumber,exppriorityVec);

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
        countt = nonbcon(model,x0_list);
        x0_list(logical(countt),:) =[];
        constantt = constantt + sum(countt);
    end
% else
%     x0_list = lhs(runmax,nParams,plb,pub,[],1e3);
% end

% optimize for starting values in RUNLIST
ML_parameters = [];
nLLVec = [];
runlist_completed = [];
for irun = 1:length(runlist)
    runlist(irun)
   
    rng(runlist(irun));
    
    x0 = x0_list(runlist(irun),:);
    fun = @(x) calc_nLL(model,x,data,fixparams,exppriorityVec);
    
    x = bads(fun,x0,lb,ub,plb,pub,@(y) nonbcon(model,y));
    fval = fun(x);
    
    x(logflag) = exp(x(logflag));
    if isempty(fixparams)
        bfp = x;
    else
        bfp = nan(1,nParams+size(fixparams,2));
        bfp(freeparamsidx) = x;
        bfp(fixparams(1,:)) = fixparams(2,:);
    end

    ML_parameters = [ML_parameters; bfp];
    nLLVec = [nLLVec fval];
    runlist_completed = [runlist_completed runlist(irun)];

end
end

