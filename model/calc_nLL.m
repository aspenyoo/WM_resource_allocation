function nLL = calc_nLL(model,Theta,data,fixparams)
%CALC_NLL calculates negative log-likelihood
% 
%   NLL = CALC_NLL(MODEL, THETA, DATA) calculates the negative
%     likelihood of DATA given MODEL and THETA.
% 
%   NLL = CALC_NLL(MODEL, THETA, DATA, FIXPARAMS) calculates the negative
%     likelihood, where THETA are free parameters and FIXPARAMS indicates
%     which parameters are fixed and what value they are fixed to. 
%
%   ================= INPUT VARIABLES ======================
% 
%   MODEL / THETA: 
%         M1, MAXIMIZING POINTS: [Jbar_total tau alpha beta]
%         M2, FLEXIBLE: [Jbar_total tau alpha beta p_high p_med]
%         M3, PROPORTIONAL: [Jbar_total tau alpha beta]
%         M4, MINIMIZING ERROR: [Jbar_total tau alpha beta gamma]
% 
%     parameter descriptions: 
%           JBAR_TOTAL: mean total amount of resources across priorities
%           TAU: second parameter of gamma noise distribution
%           ALPHA: risk preferences for post-decision wager
%           BETA: inverse noise temperature on post-decision wager
%           GAMMA: exponent for loss in Minimizing Error model
%           P_HIGH: proportion allocated to high-priority stimulus
%           P_MED: proportion allocated to medium-priority stimulus
%
%   DATA: 1 x 3 cell each containing nTrials x (1 or 2) matrix.
%     struct organization: high-, med-, low-priority stimulus trials
%     1st column: distance between target and final saccade
%     2nd column: radius of disc (ONLY IN EXP 2)
% 
%   FIXPARAMS: 2 x (number of fixed parameters) matrix, where the 1st row
%     corresponds to the linear index of each parameter to be fixed, and
%     the second row is the value of the corresponding index. default = [];
%
% 
%   ================= OUTPUT VARIABLES ================
% 
%   NLL: negative log-likelihood

% ---------------------
%      Aspen H. Yoo
%   aspen.yoo@nyu.edu
% ---------------------

if nargin < 4; fixparams = []; end

expnumber = size(data{1},2); % experiment number

% exponentiating appropriate parameters
logflag = loadconstraints(model,expnumber);

% if there are fixed parameters
if ~isempty(fixparams)
    nParams = length(Theta) + size(fixparams,2);
    nonfixedparamidx = 1:nParams;
    nonfixedparamidx(fixparams(1,:)) = [];
    
    temptheta = nan(1,nParams);
    temptheta(nonfixedparamidx) = Theta;
    temptheta(fixparams(1,:)) = fixparams(2,:);
    
    Theta = temptheta;
end

Theta(logflag) = exp(Theta(logflag));
Jbar_total = Theta(1);
tau = Theta(2);
if (expnumber == 2); alpha = Theta(3); beta = Theta(4); end


% data stuff
priorityVec = [0.6 0.3 0.1];
nPriorities = length(priorityVec);

switch model
    case 1 % optimal: maximizing points (exp 2 only)
        % calculate the proportions that maximize expected utility
        pVec = calc_optimal_pVec(Theta);
    case 2 % not optimal
        %         if sum(Theta(end-1:end))>1 % reflect over pHigh + pMed = 1 line
        %             pVec = [1-Theta(end) 1-Theta(end-1)];
        %             pVec = [pVec 1-sum(pVec)];
        %         else
        pVec = [Theta(end-1:end) 1-sum(Theta(end-1:end))];
        %         end
    case 3 % fixed
        pVec = [0.6 0.3 0.1];
    case 4 % optimal: minimizing squared error
        pVec = calc_pVec_minerror(Theta);
        %         pVec = [0.4727 0.3343 0.1930];
end

% loading vector of disc radii
[rVec] = loadvar('rVec'); % size: (1 x nrs), where nrs = 500
rVec = rVec(:); % size: (500 x 1)

if any(isinf(pVec))
    nLL = Inf;
else
    nLL = 0;
    for ipriority = 1:nPriorities
        Jbar = Jbar_total*pVec(ipriority); % Jbar for current priority condition
        
        % clear varaibles used in previous priorities (that would fuck the code
        % up if not cleared)
        clear idx1 idx2 data_r_reshaped
        
        % get subject data
        data_distance = data{ipriority}(:,1);
        nTrials = length(data_distance);
        
        % p(J|Jbar,tau)
        [JVec] = loadvar('JVec',{Jbar,tau}); % values of J
        nJs = length(JVec);
        Jpdf = gampdf(JVec,Jbar/tau,tau); % probability of that J value
        Jpdf = Jpdf./sum(Jpdf); % normalize
        %     if any(Jpdf > 1); Jpdf = Jpdf./sum(Jpdf); return; end % weird shit happens at extremely low taus, where the probability of the lowest value is >1.
        
        % p(Shat|S,J)
        Sigma = zeros(1,2,nJs*nTrials);
        Sigma(1,:,:) = sort(repmat(sqrt(1./JVec(:)),nTrials,2),'descend')'; % SDs of diagonal matrix. sigmas in descending order --> J in ascending order
        p_Shat = mvnpdf(repmat([data_distance(:) zeros(nTrials,1)],nJs,1),0,Sigma);
        p_Shat = reshape(p_Shat,nTrials,nJs)'; % nJs x nTrials
        p_Shat(p_Shat == 0) = 1e-10; % set to arbitrarily small vluae if zero
        
        % ====== Exp 2: with disc size data ======
        if (expnumber == 2)
            data_r = data{ipriority}(:,2);
            
            % p(rVec|J,beta) (a range of r to get entire probability dist)
            pdf_r = calc_pdf_r(beta, JVec, alpha); % size: (nrs x nJs)
            
            
            xdiff = diff(rVec(1:2));
            p_r = nan(nJs,nTrials);
            for iJ = 1:nJs
                pdff = pdf_r(:,iJ);
                for itrial = 1:nTrials
                    r = data_r(itrial);
                    idx1 = find((rVec- r) <= 0,1,'last');
                    idx2 = find((rVec- r) > 0,1,'first');
                    
                    p_r(iJ,itrial) = (pdff(idx2)-pdff(idx1))/xdiff.*(r-rVec(idx1)) + pdff(idx1);
                end
            end
            
            %         % calculate p(r|J,beta) (get indices of which r in rVec is closest to actual r)
            %         data_r = data_r(:)';  % size: (1 x nTrials)
            %         firstidxs = bsxfun(@(x,y) x == x(find((x-y)<=0,1,'last')),rVec,data_r); % size: (nrs x nTrials)
            %         lastidxs = bsxfun(@(x,y) x == x(find((x-y)>0,1,'first')),rVec,data_r);
            %         idx1(:,1,:) = firstidxs; % size: (nrs x 1 x nTrials)
            %         idx2(:,1,:) = lastidxs; % size: (nrs x 1 x nTrials)
            %
            %         xdiff = diff(rVec(1:2)); % size: scalar
            %         ydiff = sum(bsxfun(@times,pdf_r,idx2)) - sum(bsxfun(@times,pdf_r,idx1)); % 1 x JVec x nTrials
            %         slope = ydiff./xdiff; % 1 x JVec x nTrials
            %
            %         % linearly interpolate
            %         data_r_reshaped(1,:,:) = data_r; % 1 x 1 x nTrials
            %         p_r = squeeze(bsxfun(@times,slope,bsxfun(@minus,data_r_reshaped,sum(bsxfun(@times,rVec,idx1)))) + sum(bsxfun(@times,pdf_r,idx1))); % JVec x nTrials
        end
        
        if (expnumber == 2) % if there is disc size data
            % \int p(Shat|S,J) p(r|J) p(J) dJ
            pTrials = sum(bsxfun(@times,p_Shat.*p_r,Jpdf')); % 1 x nTrials
        else
            % \int p(Shat|S,J) p(J) dJ
            pTrials = sum(bsxfun(@times,p_Shat,Jpdf')); % 1 x nTrials
        end
        %     if ipriority == 3; Jpdf, end
        nLL = nLL - sum(log(pTrials));
        
    end
end