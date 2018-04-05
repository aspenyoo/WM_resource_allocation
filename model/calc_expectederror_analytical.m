function expectederror = calc_expectederror_analytical(Theta,allocatedpriorityVec,exppriorityVec)
%CALC_EXPECTEDERROR_ANALYTICAL(THETA,ALLOCATEDPRIORITYVEC) calculates the
% expected cost (euclidean error ^ psi) for a given parameter vector THETA
% and resource allocation across consitions ALLOCATEDPRIORITYVEC.
% 
% ===== INPUT VARIABLES =====
% THETA: Jbar_total, tau, [alpha, beta,] gamma. (alpha and beta are only in
%   experiment 2)
% ALLOCATEDPRIORITYVEC: row vector of allocated priority.
%   sum(allocatedpriorityVec) = 1
% EXPPRIORITYVEC: row vector of experimental priority.
%   sum(exppriorityVec) = 1

%if nargin < 3; exppriorityVec = [0.6 0.3 0.1]; end

% getting parameters
Jbar_total = Theta(1);
tau = Theta(2);
gamma = Theta(end);

nPriorities = length(allocatedpriorityVec);

% E[d^gamma] = \int d^blah \int p(d|J)p(J) dd dJ
% E[d^gamma] = \int d^blah \int rayleigh(d|1/J) gamma(J|Jbar,tau) dd dJ

expectederror = 0;
for ipriority = 1:nPriorities
    Jbar = Jbar_total*allocatedpriorityVec(ipriority);

    k = Jbar/tau;
    
    if any([gamma/2+1 k-(gamma/2) k] <= 0) % this constraint must be met for expected error to be analytical
        % ASPEN, why not do something not analytical here? 
        expectederror = Inf;
        return
    else
        bleh = exp(gammaln(gamma/2 + 1) + gammaln(k-(gamma/2)) - gammaln(k)) .* (2/tau)^(gamma/2);
    end
    
    expectederror = expectederror + exppriorityVec(ipriority).*bleh;
    
end

