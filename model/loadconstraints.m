function [logflag, lb, ub, plb, pub, nonbcon] = loadconstraints(model,exppriorityVec,is_wagerdata)
%LOADCONSTRAINTS loads optimization constraints for each model
% 
% LOADCONSTRAINTS(MODEL,EXPPRIORITYVEC) loads optimization constraint for a
% particular model
%
% LOADCONSTRAINTS(MODEL,EXPPRIORITYVEC,IS_WAGERDATA) loads optimization
% constraint for wager data, if IS_WAGERDATA = 1
% 
% ===== INPUT VARIABLES =====
% 
% MODEL: 'max_points', 'flexible', 'proportional', or 'min_error'
% 
% EXPPRIORITYVEC: row vector of experimental priority.
%   sum(exppriorityVec) = 1
%
% IS_WAGERDATA: 0: just error data. 1: error and wager data

if nargin < 3; is_wagerdata = 0; end

% change values if fitting all conditions jointly
if strcmp(model(1:3),'all'); 
    model = model(5:end); 
    allflag = 1;
else
    allflag = 0;
end

nPriorities = length(exppriorityVec); 

% lower and upper bounds, logflags
lb = [1e-5 1e-3]; % Jbar_total, tau
ub = [50 10];
plb = [0.5 0.01];
pub = [10 1];
logflag = [1 1];

if (allflag) % if fitting all TMS conditions jointly
    lb = lb([1 1 2]); % Jbar_total, Jbar_total, tau
    ub = ub([1 1 2]);
    plb = plb([1 1 2]);
    pub = pub([1 1 2]);
    logflag = logflag([1 1 2]);
end

if (is_wagerdata) % alpha beta
    lb = [lb 1e-5 1e-5];
    ub = [ub 5 5];
    plb = [plb 0.7 0.5];
    pub = [pub 1.3 1.5];
    logflag = [logflag 0 0];
end

switch model
    case {'flexible','stronghyp'} % define p_high p_med for Flexible model
        lb = [lb 1e-10.*ones(1,nPriorities-1)];
        ub = [ub ones(1,nPriorities-1)-eps];
        plb = [plb max([1e-10.*ones(1,nPriorities-1); exppriorityVec(1:end-1).*0.5])];
        pub = [pub min([ones(1,nPriorities-1)-eps; exppriorityVec(1:end-1).*1.5])];
        logflag = [logflag zeros(1,nPriorities-1)];
        
        if (allflag) % if fitting all TMS conditions jointly
            lb = [lb lb(end).*ones(1,nPriorities-1)]; 
            ub = [ub ub(end).*ones(1,nPriorities-1)];
            plb = [plb plb(end).*ones(1,nPriorities-1)];
            pub = [pub pub(end).*ones(1,nPriorities-1)];
            logflag = [logflag logflag(end).*ones(1,nPriorities-1)];
        end
        
    case 'min_error' % gamma for Minimizing Error model
        lb = [lb 1e-10];
        ub = [ub 10];
        plb = [plb 1e-3];
        pub = [pub 1];
        logflag = [logflag 1];
        nonbcon = @model4nonbcon; % violates if Jbar/tau - gamma/2 <= 0 
    otherwise
        nonbcon = [];
end

switch model
    case 'stronghyp' % adds another jbar and p's to fit
        lb = [lb 1e-5 1e-10.*ones(1,nPriorities-1)];
        ub = [ub 50 ones(1,nPriorities-1)-eps];
        plb = [plb 0.5 max([1e-10.*ones(1,nPriorities-1); exppriorityVec(1:end-1).*0.5])];
        pub = [pub 10 min([ones(1,nPriorities-1)-eps; exppriorityVec(1:end-1).*1.5])];
        logflag = [logflag 1 zeros(1,nPriorities-1)];
end

logflag = logical(logflag); 
lb(logflag) = log(lb(logflag)); 
ub(logflag) = log(ub(logflag)); 
plb(logflag) = log(plb(logflag)); 
pub(logflag) = log(pub(logflag)); 

