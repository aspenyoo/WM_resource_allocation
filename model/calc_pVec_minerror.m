function [pVec, fval] = calc_pVec_minerror(Theta)
% PVEC = calc_pVec_optimalerror(THETA)
% 
% calculates the proportion allocated to each priority condition that
% minimizes the squared error of target and response. 
Theta

% function for expected squared error
calc_E_error = @(x) calc_expectederror_analytical(Theta,x);

% parameters for optimization
Aeq = [1 1 1];
beq = 1;
A = diag(ones(1,3).*(-Theta(1)/Theta(2)));
b = -Theta(end)/2.*ones(3,1);
nonlcon = deal([]);
options = optimset('Display','none');
lb = [1e-3 1e-3 1e-3];
ub = [1 1 1];
nStartVals = 10;  % tried with different parameters and lowest value showed up 3,5,7,8,10,10 of 10. 

x0 = rand(nStartVals,2);
x0 = [x0 1-sum(x0,2)];

while any((x0(:,3) < 0) | sum(A*x0' > -(Theta(end)/2))') % make sure it satisfies Aeq beq constraints
    idx = (x0(:,3) < 0) | sum(A*x0' > -(Theta(end)/2))';
    x0(idx,1:2) = rand(sum(idx),2);
    x0(idx,3) = 1 - sum(x0(idx,1:2),2);
end

% optimizing
pVec = nan(nStartVals,3);
nEU = nan(1,nStartVals);
for istartval = 1:nStartVals
    [pVec(istartval,:), nEU(istartval)] = fmincon(calc_E_error,x0(istartval,:),A,b,Aeq,beq,lb,ub,nonlcon,options);
end

pVec(nEU < 0,:) = [];
nEU(nEU < 0) = [];

fval = min(nEU);
pVec = pVec(nEU == fval,:);

if isempty(fval)
    fval = Inf;
    pVec = Inf*ones(1,3);
else
    pVec = pVec(1,:);  % in case multiple entries have the nEU == min(nEU)
end

