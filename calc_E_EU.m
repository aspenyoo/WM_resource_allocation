function E_EU = calc_E_EU(Theta)
% CALC_E_EU calculates expected value of the expected utility for a given 
% Theta = [Jbar, tau, beta]. 
%
% aspen yoo -- 02.08.2016

% % input parameters
% Theta = [15 1 1];

Jbar = Theta(1);
tau = Theta(2);
beta = Theta(3);

% loading in variables
[JVec,rVec] = loadvar('JVec','rVec');

% p(r|theta)
p_r = calc_pdf_r(beta,JVec);

% EU(J) = \int EU(r,J) p(r) dr
EU_J = sum(calc_EU(rVec,JVec).*p_r);

% for gamma distribution
Jpdf = gampdf(JVec,Jbar/tau,tau);
Jpdf = Jpdf./sum(Jpdf);

% EU(Jbar) = \int EU(J) p(J) dJ
E_EU = Jpdf*EU_J';
