function [data] = simulate_data(model,expnumber,Theta,nTrials,expPriorityVec)
% SIMULATE_DATA simulates data for the optimal observer in the [0.6 0.3
%   0.1] priority task.
%
% ===================== INPUT VARIABLES ========================
% THETA: 1 x N vector of parameters: Jbar_total, tau[, alpha, beta, gamma]
% EXPNUMBER: 1 (no disc). 2 (with disc).
% NTRIALS: 1 x 3 vector of number of trials for each priority value [0.6
%   0.3 0.1]
%
% ===================== OUTPUT VARIABLE =======================
% DATA: 1 x 3 cell (for each priority), each containing a
%   nTrials(ipriority) x 2 matrix where column 1 is the error distance (dva)
%   and the second column is the size of the disc (dva)
%
%
% -----------------------
%      Aspen H. Yoo
%   aspen.yoo@nyu.edu

% if nargin < 4; nTrials = [250 120 70]; end % mean number of trials across actual participants
% if nargin < 5; expPriorityVec = [0.6 0.3 0.1]; end

nPriorities = length(expPriorityVec);
nSubj = size(Theta,1);

data = cell(1,nSubj);
for isubj = 1:nSubj
    theta = Theta(isubj,:);
    
    Jbar_total = theta(1);
    tau = theta(2);
    if (expnumber == 2)
        alpha = theta(3);
        beta = theta(4);
        rVec = loadvar('rVec');
    end
    
    switch model
        case 'max_points'
            pVec = calc_optimal_pVec(theta,expPriorityVec);
        case 'flexible'
            pp = Theta(end-(nPriorities-2):end);
            pVec = [pp 1-sum(pp)];
        case 'proportional'
            pVec = expPriorityVec;
        case 'min_error' 
            pVec = calc_pVec_minerror(theta,expPriorityVec);
    end
    
    % make data
    data{isubj} = cell(1,nPriorities);
    for ipriority = 1:nPriorities
        clear Sigma
        
        Jbar = Jbar_total*pVec(ipriority); % mean precision
        ntrials = nTrials(ipriority); % number of trials
        data{isubj}{ipriority} = nan(ntrials,expnumber);
        
        JVec = gamrnd(Jbar/tau,tau,[1 ntrials]); % precision on each trial
        JVec(JVec < 1e-10) = 1e-10;
        
        % generating Shat
        Sigma(1,:,:) = repmat(sqrt(1./JVec(:)),1,2)'; % get diagonals of covariance matrix
        errros = mvnrnd([0 0],Sigma); % generating data
        data{isubj}{ipriority}(:,1) = sqrt(sum(errros.^2,2));
        
        if (expnumber == 2)
            % generating disc size
            pdf_r = calc_pdf_r(beta, JVec, alpha); % length(rVec) x length(JVec)
            cdf_r = cumsum(pdf_r);
            cdf_r = cdf_r./cdf_r(end); % normalize
            samples = num2cell(rand(1,ntrials));
            cdf_r = num2cell(cdf_r,1);
            
            idxs = cell2mat(cellfun(@(x,y) find(x>=y,1,'first'),cdf_r,samples,'UniformOutput',false));
            data{isubj}{ipriority}(:,2) = rVec(idxs);
        end
    end
end


end