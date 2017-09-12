% function calc_expectederror(Theta,priorityVec)

% getting parameters
Jbar_total = Theta(1);
tau = Theta(2);
if (expnumber == 2); alpha = Theta(3); beta = Theta(4); end
if (model == 4); blah = Theta(end); end

priorityVec = [0.1 0.3 0.6];
nPriorities = length(priorityVec);

for ipriority = 1:nPriorities
    Jbar = Jbar_total*priorityVec(ipriority);
    
    % get x-axis of J and distance d 
    [JVec,dVec] = loadvar('JVec',{Jbar,tau},'rVec');
    dVec = dVec';
    
    % get Jpdf
    nJs = length(JVec);
    Jpdf = gampdf(JVec,Jbar/tau,tau); % probability of that J value
    Jpdf = Jpdf./sum(Jpdf);
    
    % 
    
end