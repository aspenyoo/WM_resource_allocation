function fit_parameters(subjnum,nStartVals)
if nargin < 2; nStartVals = 1; end

load('cleandata.mat')
subjdata = data{subjnum};

filepath = '/home/ay963/spatialWM/fits/';
filename = [filepath 'fits_optimal_subj' num2str(subjnum) '.mat'];

lb = [1e-5 1e-5 1e-5]; % Jbar_total, tau, beta, lapse (ASPEN FIGURE OUT LAPSE STUFF)
ub = [50 50 10]; % ASPEN: refine
plb = [1e-5 1e-3 1e-5];
pub = [30 10 5]; % ASPEN: refine
logflag = logical([1 1 0]);
lb(logflag) = log(lb(logflag));
ub(logflag) = log(ub(logflag));
plb(logflag) = log(plb(logflag));
pub(logflag) = log(pub(logflag));

for istartvals = 1:nStartVals;
    try load(filename); catch; ML_parameters = []; nLLVec = []; end
    startingval = plb + rand(1,3).*(pub - plb);
    objfunc = @(x) calc_nLL(subjdata,x);
    [x,fval] = bps(objfunc,startingval,lb,ub,plb,pub);
    x(logflag) = exp(x(logflag));
    ML_parameters = [ML_parameters; x];
    nLLVec = [nLLVec fval];
    save(filename,'ML_parameters','nLLVec')
end