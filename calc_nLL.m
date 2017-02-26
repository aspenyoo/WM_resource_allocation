function nLL = calc_nLL(model,Theta,data)
% CALC_NLL(JBAR_TOTAL,TAU,BETA)
%
% CALC_NLL: calculates negative log likelihood of parameter combination for
% the optimal model
% 
% ================= INPUT VARIABLES ======================
% MODEL: 1 (optimal) or 2 (not optimal)
% 
% THETA = vector of parameters. For model 1, [ Jbar_total tau beta]. 
%   For model 2, [Jbar_total tau beta p_high p_med]
%       JBAR_TOTAL: total amount of resources across priorities
%       TAU: second parameter of gamma noise distribution
% 
% DATA: 1 x 3 cell each containing nTrials x 2 matrix. 
%   struct organization: high, med, low priority trials
%   1st column: distance between target and final saccade
%   2nd column: radius of disc
% 
% ================= OUTPUT VARIABLES ================
% NLL: -L(Theta|data,model)
%
% -----------------------
%      Aspen H. Yoo
%   aspen.yoo@nyu.edu


% Theta = [5 1 1];
switch model
    case 1  % optimal model
        logflag = logical([1 1 0]);
    case 2 % not optimal model
        logflag = logical([1 1 0 0 0]);
        
        % if the sum of the proportions allocated to high and medium
        % priorities are greater than 1
        if sum(Theta(4:5)) > 1
            nLL = Inf;
            return;
        end
end
Theta(logflag) = exp(Theta(logflag));
Jbar_total = Theta(1);
tau = Theta(2);
beta = Theta(3);
% lapse = Theta(4);


% data stuff
priorityVec = [0.6 0.3 0.1];
nPriorities = length(priorityVec);

switch model
    case 1 % optimal
        
        % calculate the optimal proportions given the parameters
        calc_ntotalEU = @(x) -(0.6*calc_E_EU([Jbar_total*x(1),tau,beta]) ...
            + 0.3*calc_E_EU([Jbar_total*x(2),tau,beta])...
            + 0.1*calc_E_EU([Jbar_total*x(3),tau,beta]));
        
        Aeq = [1 1 1];
        beq = 1;
        [A,b,nonlcon] = deal([]);
        options = optimset('Display','none');
        lb = [1e-5 1e-5 1e-5];
        ub = [1 1 1];
        nStartVals = 5; % tried with different parameters and lowest value showed up 3,5,7,8,10,10 of 10.
        pVec = nan(nStartVals,3);
        nEU = nan(1,nStartVals);
        for istartval = 1:nStartVals
            [pVec(istartval,:), nEU(istartval)] = fmincon(calc_ntotalEU,rand(1,3),A,b,Aeq,beq,lb,ub,nonlcon,options);
        end
        pVec = pVec(nEU == min(nEU),:);
        pVec = pVec(1,:); % in case multiple entries have the nEU == min(nEU)
    
    case 2 % not optimal
        pVec = [Theta(4:5) 1-sum(Theta(4:5))];
end

% loading rVec;
[rVec] = loadvar('rVec');
rVec = rVec(:); % vertical

nLL = 0;
for ipriority = 1:nPriorities
    ipriority
    Jbar = Jbar_total*pVec(ipriority); % Jbar for current trial
    
    % get data
    data_distance = data{ipriority}(:,1);
    data_r = data{ipriority}(:,2);
    
    % p(J|Jbar,tau)
    [JVec] = loadvar({'JVec',Jbar,tau});
    nJs = length(JVec);
    Jpdf = gampdf(JVec,Jbar/tau,tau);
    Jpdf = Jpdf./qtrapz(Jpdf); % normalize
    
    % p(Shat|S,J)
    nTrials = length(data_distance);
    Sigma = zeros(1,2,nJs*nTrials);
    Sigma(1,:,:) = sort(repmat(sqrt(1./JVec(:)),nTrials,2),'descend')'; % in descending order to keep J ascending
    p_Shat = mvnpdf(repmat([data_distance(:) zeros(nTrials,1)],nJs,1),0,Sigma);
    p_Shat = reshape(p_Shat,nTrials,nJs);
    
    % p(Shat|S) = \int p(Shat|S,J) p(J) dJ
    p_Shat = qtrapz(bsxfun(@times,p_Shat,Jpdf),2);
    
    % get pdf of p(r|Jbar,tau,beta)
    pdf_r = calc_pdf_r(beta, JVec);
    pdf_r = bsxfunandsum(@times,pdf_r,Jpdf,2);
    
    % p(r): probability of responding r given the parameters
    data_r = data_r(:)';  % make sure it is horizontal vector
    size(data_r)
    firstidxs = bsxfun(@(x,y) x == x(find((x-y)<=0,1,'last')),rVec,data_r);
    lastidxs = bsxfun(@(x,y) x == x(find((x-y)>=0,1,'first')),rVec,data_r);
    
    % linearly interpolate
    slope = (sum(bsxfun(@times,pdf_r,lastidxs)) - sum(bsxfun(@times,pdf_r,firstidxs)))./diff(rVec(1:2));
    p_r = slope.*(data_r - sum(bsxfun(@times,rVec,firstidxs))) + sum(bsxfun(@times,pdf_r,firstidxs));
    
    % nLL
    nLL = nLL -sum(log(p_Shat)) -sum(log(p_r));
    %     nLL = nLL -sum(log((1-lapse).*p_Shat + lapse.*p_Shat_lapse)) -sum(log(p_r));
end