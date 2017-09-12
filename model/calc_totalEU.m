function totalEU = calc_totalEU(Theta, model)
% CALC_TOTALEU calculates expected value of the expected utility across the
% entire experiment
%
% aspen yoo -- 08.21.2017

Jbar_total = Theta(1);
tau = Theta(2);
alpha = Theta(3);
beta = Theta(4);

% get pVec
switch model
    case 1 
        pVec = calc_optimal_pVec(Theta);
    case 2
        pVec = [Theta(5:6) 1-sum(Theta(5:6))];
    case 3
        pVec = [0.6 0.3 0.1];
    case 4
        pVec = [0.4727 0.3343 0.1930];
end

totalEU = 0.6*calc_E_EU([Jbar_total*pVec(1),tau,alpha,beta]) ...
    + 0.3*calc_E_EU([Jbar_total*pVec(2),tau,alpha,beta])...
    + 0.1*calc_E_EU([Jbar_total*pVec(3),tau,alpha,beta]);