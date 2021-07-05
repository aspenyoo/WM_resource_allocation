function nLL = calc_nLL_allcond(model,Theta,data,exppriorityVec,fixparams)

% hypothesis is that tms to sPCS decreases priority effect
%                        to IPS2 lowers overall Jbartotal

isviolated = check_nonbcon_allcond(model,Theta,exppriorityVec);
if (isviolated)
    nLL = Inf;
else
    % Theta = [Jbar_total_noTMS Jbar_total_ips2 tau phigh_noTMS phigh_spcs]
    
    model = model(5:end);
    fp = [];
    
    % noTMS condition
    d = data.noTMS;
    x = Theta([1 3 4]);
    if ~isempty(fixparams); fp = fixparams([1 3 4]); end
    nLL = calc_nLL(model,x,d,exppriorityVec,fp);
    
    % l IPS2
    d = data.l_ips2;
    x = Theta([2 3 4]);
    if ~isempty(fixparams); fp = fixparams([2 3 4]); end
    nll = calc_nLL(model,x,d,exppriorityVec,fp);
    nLL = nLL+nll;s
    
    % l SPCS 
    d = data.l_spcs;
    x = Theta([1 3 5]);
    if ~isempty(fixparams); fp = fixparams([1 3 5]); end
    nll = calc_nLL(model,x,d,exppriorityVec,fp);
    nLL = nLL+nll;
end
