function E_EU = calc_E_EU(Theta)
%CALC_E_EU calculates expected value of the expected utility for a given 
% Theta = [Jbar, tau, alpha, beta] for Maximizing Points model

% ---------------------
%      Aspen H. Yoo
%   aspen.yoo@nyu.edu
% ---------------------

Jbar = Theta(1);   % mean of memory precision gamma distribution
tau = Theta(2);    % scale parameter of memeory precision gamma distribution
alpha = Theta(3);  % risk preferences for post-decision wager
beta = Theta(4);   % inverse noise temperature on post-decision wager

% loading in variables
[JVec,rVec] = loadvar('JVec',{Jbar,tau},'rVec');

% p(r|beta,J)
p_r = calc_pdf_r(beta,JVec,alpha);

% EU(J) = \int EU(r,J) p(r) dr
EU_J = sum(calc_EU(rVec,JVec,alpha).*p_r);

% p(J). J ~ Gamma(Jbar, tau)
Jpdf = gampdf(JVec,Jbar/tau,tau);
Jpdf = Jpdf./sum(Jpdf);

% EU(Jbar) = \int EU(J) p(J) dJ
E_EU = sum(Jpdf.*EU_J);
