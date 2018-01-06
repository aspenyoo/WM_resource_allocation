function expectederror = calc_expectederror_analytical(Theta,allocatedpriorityVec)
% CALC_EXPECTEDERROR_ANALYTICAL(THETA,ALLOCATEDPRIORITYVEC) calculates the
% expected cost (euclidean error ^ psi) for a given parameter vector THETA
% and resource allocation across consitions ALLOCATEDPRIORITYVEC.
% 
% ===== INPUT VARIABLES =====
% THETA: Jbar_total, tau, (alpha, beta,) psi (alpha and beta are only in
%   experiment 2)
% ALLOCATEDPRIORITYVEC: 1 x 3 vector of allocated priority.
%   sum(allocatedpriorityVec) = 1

% getting parameters
Jbar_total = Theta(1);
tau = Theta(2);
gamma = Theta(end);

priorityVec = [0.6 0.3 0.1];
nPriorities = length(allocatedpriorityVec);

% E[d^gamma] = \int d^blah \int p(d|J)p(J) dd dJ
% E[d^gamma] = \int d^blah \int rayleigh(d|1/J) gamma(J|Jbar,tau) dd dJ

expectederror = 0;
for ipriority = 1:nPriorities
    Jbar = Jbar_total*allocatedpriorityVec(ipriority);

    k = Jbar/tau;
    
    if any([gamma/2+1 k-(gamma/2) k] <=0)
        expectederror = Inf;
        return
    else
        bleh = exp(gammaln(gamma/2 + 1) + gammaln(k-(gamma/2)) - gammaln(k)).* (2/tau)^(gamma/2);
    end
    
    expectederror = expectederror + priorityVec(ipriority).*bleh;
    
end