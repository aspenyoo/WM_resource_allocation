function nLL = calc_nLL(model,Theta,data,fixparams,exppriorityVec)
%CALC_NLL calculates negative log-likelihood of parameters given data and
%model.
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
%         'max_points', MAXIMIZING POINTS / [Jbar_total tau alpha beta]
%         'flexible', FLEXIBLE / [Jbar_total tau alpha beta p_high p_med]
%         'proportional', PROPORTIONAL / [Jbar_total tau alpha beta]
%         'min_error', MINIMIZING ERROR / [Jbar_total tau alpha beta gamma]
% 
%       parameter descriptions: 
%           JBAR_TOTAL: mean total amount of resources across priorities
%           TAU: second parameter of gamma noise distribution
%           ALPHA: risk preferences for post-decision wager
%           BETA: inverse noise temperature on post-decision wager
%           GAMMA: exponent for loss in Minimizing Error model
%           P_HIGH: proportion allocated to high-priority stimulus
%           P_MED: proportion allocated to medium-priority stimulus
%
%   DATA: cell of length nPriorities. the ith element of DATA should be
%     all data corresponding to EXPPRIORITYVEC(i) condition. the first
%     column should contain the magnitude of errors and the second column,
%     if available, should contain the corresponding circle wager radius
%     size. 
% 
%   FIXPARAMS: (optional). 2 x (number of fixed parameters) matrix. fixed 
%     parameters, such that the first row corresponds to the index and 
%     second row corresponds to the value of the fixed parameter. 
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
logflag = loadconstraints(model,exppriorityVec,expnumber-1);

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

% obtain vector of resource allocated
nPriorities = length(exppriorityVec);
switch model
    case 'max_points'   % maximizing points (exp 2 only)
        pVec = calc_pVec_maxpoints(Theta,exppriorityVec);
    case 'flexible'     % flexible
        pp = Theta(end-(nPriorities-2):end);
        pVec = [pp 1-sum(pp)];
    case 'proportional' % proportional
        pVec = exppriorityVec;
    case 'min_error'    % minimizing error^\gamma
        pVec = calc_pVec_minerror(Theta,exppriorityVec);
end

% loading vector of circle wager radii
[rVec] = loadvar('rVec'); % size: (1 x 500)
rVec = rVec(:); % size: (500 x 1)

if any(isinf(pVec))
    nLL = Inf;
else
    nLL = 0;
    for ipriority = 1:nPriorities
        Jbar = Jbar_total*pVec(ipriority); % Jbar for current priority condition
        
        % clear varaibles used in previous priorities (necessary for code to run)
        clear idx1 idx2 data_r_reshaped
        
        % get subject data
        data_distance = data{ipriority}(:,1);
        nTrials = length(data_distance);
        
        % p(J|Jbar,tau)
        [JVec] = loadvar('JVec',{Jbar,tau}); % values of J
        nJs = length(JVec);
        Jpdf = gampdf(JVec,Jbar/tau,tau); % probability of that J value
        Jpdf = Jpdf./sum(Jpdf); % normalize
        
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
        end
        
        if (expnumber == 2) % if there is wager data
            % \int p(Shat|S,J) p(r|J) p(J) dJ
            pTrials = sum(bsxfun(@times,p_Shat.*p_r,Jpdf')); % 1 x nTrials
        else
            % \int p(Shat|S,J) p(J) dJ
            pTrials = sum(bsxfun(@times,p_Shat,Jpdf')); % 1 x nTrials
        end
        
        nLL = nLL - sum(log(pTrials));
        
    end
end