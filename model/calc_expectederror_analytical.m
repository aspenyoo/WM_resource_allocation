function expectederror = calc_expectederror_analytical(Theta,allocatedpriorityVec)

% allocatedpriorityVec

% getting parameters
Jbar_total = Theta(1);
tau = Theta(2);
psi = Theta(end);

priorityVec = [0.6 0.3 0.1];
nPriorities = length(allocatedpriorityVec);

% E[d^gamma] = \int d^blah \int p(d|J)p(J) dd dJ
% E[d^gamma] = \int d^blah \int rayleigh(d|1/J) gamma(J|Jbar,tau) dd dJ

expectederror = 0;
for ipriority = 1:nPriorities
    Jbar = Jbar_total*allocatedpriorityVec(ipriority);

    k = Jbar/tau;
    
    bleh = exp(gammaln(psi/2 + 1) + gammaln(k-(psi/2)) - gammaln(k)).* (2/tau)^(psi/2);
    
    expectederror = expectederror + priorityVec(ipriority).*bleh;
    
end