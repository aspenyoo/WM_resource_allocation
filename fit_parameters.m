function fit_parameters(subjnum,nStartVals)
if nargin < 2; nStartVals = 1; end

load('cleandata.mat')
subjdata = data{subjnum};

filepath = '/home/ay963/spatialWM/fits/';
filename = [filepath 'fits_optimal_subj' num2str(subjnum) '.mat'];

lb = [1e-5 1e-5 1e-5]; % Jbar_total, tau, beta, lapse (ASPEN FIGURE OUT LAPSE STUFF)
ub = [50 50 5]; % ASPEN: refine
plb = [1e-5 1e-3 1e-5];
pub = [30 10 2]; % ASPEN: refine
logflag = logical([1 1 0]);
lb(logflag) = log(lb(logflag));
ub(logflag) = log(ub(logflag));
plb(logflag) = log(plb(logflag));
pub(logflag) = log(pub(logflag));

for istartvals = 1:nStartVals;
    try load(filename); catch; ML_parameters = []; nLLVec = []; end
    x0 = plb + rand(1,3).*(pub - plb);
    [x,fval] = bps(@(x) calc_nLL(x,subjdata),x0,lb,ub,plb,pub);
    x(logflag) = exp(x(logflag));
    ML_parameters = [ML_parameters; x];
    nLLVec = [nLLVec fval];
    save(filename,'ML_parameters','nLLVec')
end