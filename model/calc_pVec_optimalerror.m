function pVec = calc_pVec_optimalerror(Theta)
% PVEC = calc_pVec_optimalerror(THETA)
% 
% calculates the proportion allocated to each priority condition that
% minimizes the squared error of target and response. 

Jbar_total = Theta(1);
tau = Theta(2);
alpha = Theta(3);
beta = Theta(4);

% data stuff
priorityVec = [0.6 0.3 0.1];
nPriorities = length(priorityVec);

% equation for expected squared error
calc_E_squarederror = @(x) 0.6/x(1) + 0.3/x(2) + 0.1/x(3);


Aeq = [1 1 1];
beq = 1;
[A,b,nonlcon] = deal([]);
options = optimset('Display','none');
lb = [1e-5 1e-5 1e-5];
ub = [1 1 1];
nStartVals = 10; % tried with different parameters and lowest value showed up 3,5,7,8,10,10 of 10. 
pVec = nan(nStartVals,3);
nEU = nan(1,nStartVals);
for istartval = 1:nStartVals
    [pVec(istartval,:), nEU(istartval)] = fmincon(calc_E_squarederror,rand(1,3),A,b,Aeq,beq,lb,ub,nonlcon,options);
end
pVec = pVec(nEU == min(nEU),:);
pVec = pVec(1,:); % in case multiple entries have the nEU == min(nEU)
