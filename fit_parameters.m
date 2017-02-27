function fit_parameters(model,subjnum,nStartVals)
if nargin < 3; nStartVals = 1; end
%
%
% ================= INPUT VARIABLES ==================
% MODEL: 1 (optimal priority placement) or 2 (not optimal)
% SUBJNUM: subject number. 1 - 11
% NSTARTVALS: optimization starting vals
% 
% -----------------------
%      Aspen H. Yoo
%   aspen.yoo@nyu.edu


load('cleandata.mat')
subjdata = data{subjnum};

% filepath = 'fits/';
filepath = '/home/ay963/spatialWM/fits/';
filename = [filepath 'fits_model' num2str(model) '_subj' num2str(subjnum) '.mat'];

lb = [1e-5 1e-5 1e-5]; % Jbar_total, tau, beta, lapse (ASPEN FIGURE OUT LAPSE STUFF)
ub = [50 50 5]; % ASPEN: refine
plb = [1e-5 1e-3 1e-5];
pub = [30 10 2]; % ASPEN: refine
logflag = logical([1 1 0]);
if model == 2
    lb = [lb 0 0];
    ub = [ub 1 1];
    plb = [plb 0 0];
    pub = [pub 1 1];
    logflag = logical([logflag 0 0]);
end
nParams = length(logflag);
lb(logflag) = log(lb(logflag));
ub(logflag) = log(ub(logflag));
plb(logflag) = log(plb(logflag));
pub(logflag) = log(pub(logflag));

for istartvals = 1:nStartVals
    try load(filename); catch; ML_parameters = []; nLLVec = []; end
    x0 = plb + rand(1,nParams).*(pub - plb);
    [x,fval] = bps(@(x) calc_nLL(model,x,subjdata),x0,lb,ub,plb,pub);
    x(logflag) = exp(x(logflag));
    ML_parameters = [ML_parameters; x];
    nLLVec = [nLLVec fval];
    save(filename,'ML_parameters','nLLVec')
end