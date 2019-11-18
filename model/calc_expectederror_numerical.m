function expectederror = calc_expectederror_numerical(Theta,allocatedpriorityVec,exppriorityVec)

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
    [JVec,dVec] = loadvar('JVec',{Jbar,tau},'rVec');
    JVec = JVec';
    
    % get Jpdf: p(J|Jbar,tau). 500 x 1 vector
    Jpdf = gampdf(JVec,Jbar/tau,tau); % probability of that J value
    Jpdf = Jpdf./sum(Jpdf);
    
    % get dpdf given J: p(d|1/J). will be nJs (500) x nds (500) vector
    d_given_J = bsxfun(@(x,y) x.*y.*exp(-x^2.*y/2),dVec,JVec);
    
    % get dpdf (marginalize over J)
    dpdf = bsxfunandsum(@times,d_given_J,Jpdf);
    dpdf = dpdf./sum(dpdf);
    
    % \int d^blah p(d) dd
    expectederror = expectederror + exppriorityVec(ipriority).*sum((dVec.^gamma) .* dpdf);
end