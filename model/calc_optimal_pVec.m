function pVec = calc_optimal_pVec(Theta)

Jbar_total = Theta(1);
tau = Theta(2);
alpha = Theta(3);
beta = Theta(4);

% data stuff
priorityVec = [0.6 0.3 0.1];
nPriorities = length(priorityVec);

% calculate the optimal proportions given the parameters
calc_ntotalEU = @(x) -(0.6*calc_E_EU([Jbar_total*x(1),tau,alpha,beta]) ...
    + 0.3*calc_E_EU([Jbar_total*x(2),tau,alpha,beta])...
    + 0.1*calc_E_EU([Jbar_total*x(3),tau,alpha,beta]));

Aeq = [1 1 1];
beq = 1;
[A,b,nonlcon] = deal([]);
options = optimset('Display','none');
lb = [1e-5 1e-5 1e-5];
ub = [1 1 1];
nStartVals = 5; % tried with different parameters and lowest value showed up 3,5,7,8,10,10 of 10. 
pVec = nan(nStartVals,3);
nEU = nan(1,nStartVals);
for istartval = 1:nStartVals
    [pVec(istartval,:), nEU(istartval)] = fmincon(calc_ntotalEU,rand(1,3),A,b,Aeq,beq,lb,ub,nonlcon,options);
end
pVec = pVec(nEU == min(nEU),:);
pVec = pVec(1,:); % in case multiple entries have the nEU == min(nEU)
