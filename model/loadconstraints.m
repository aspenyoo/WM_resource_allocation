function [logflag, lb, ub, plb, pub, nonbcon] = loadconstraints(model,expnumber)

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
switch model
    case 2 % define p_high p_med for Flexible model
        lb = [lb 1e-10 1e-10];
        ub = [ub 1 1];
        plb = [plb 0.3 1e-10];
        pub = [pub 0.7 0.3];
        logflag = [logflag 0 0];
    case 4 % gamma for Minimizing Error model
        lb = [lb 1e-10];
        ub = [ub 10];
        plb = [plb 1e-3];
        pub = [pub 1];
        logflag = [logflag 1];
        nonbcon = @model4nonbcon; % violates if Jbar/tau - gamma/2 <= 0 
    otherwise
        nonbcon = [];
end
logflag = logical(logflag); 
lb(logflag) = log(lb(logflag)); 
ub(logflag) = log(ub(logflag)); 
plb(logflag) = log(plb(logflag)); 
pub(logflag) = log(pub(logflag)); 


