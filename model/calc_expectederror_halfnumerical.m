function expectederror = calc_expectederror_halfnumerical(Theta,allocatedpriorityVec,exppriorityVec)
%CALC_EXPECTEDERROR_NUMERICAL numerically computes the expected error of a
%given alllocatedpriorityVec
%
% 10/25/2018: the current implementation here is COMPLETELY INCORRECT!!

% -----------------------
%      Aspen H. Yoo
%   aspen.yoo@nyu.edu
% -----------------------

% getting parameters
Jbar_total = Theta(1);
tau = Theta(2);
gamma = Theta(end);

nPriorities = length(exppriorityVec);

% E[d^gamma] = \int d^blah \int p(d|J)p(J) dd dJ
% E[d^gamma] = \int d^blah \int rayleigh(d|1/J) gamma(J|Jbar,tau) dd dJ

expectederror = 0;
for ipriority = 1:nPriorities
    Jbar = Jbar_total*allocatedpriorityVec(ipriority);
    
    % get x-axis of J and distance d 
    [JVec] = loadvar('JVec',{Jbar,tau});
    dJ = diff(JVec(1:2));
    JVec = JVec';
    
    % get value of scalar
    k = Jbar/tau;
    scalar = exp(gammaln(gamma/2+1) - gammaln(k))*(2^(gamma/2)/tau^k);
    blah = exppriorityVec(ipriority).*sum(JVec.^(k-gamma/2-1).*exp(-k))*dJ;
    
    % \int d^blah p(d) dd
    expectederror = expectederror + blah;
    expectederror
end