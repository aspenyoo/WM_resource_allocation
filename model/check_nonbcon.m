function isviolated = check_nonbcon(model,x,exppriorityVec)
%CHECK_NONBCON checks for non-bound constraints
%
%   CHECK_NONBCON checks if X violates any of the non-bound constraints for
%     MODEL. 
%
%   ========== INPUT VARIABLES ==========
% 
%   MODEL: 'proportional', 'min_error', 'flexible', 'max_points'
%
%   X: parameter vector

isviolated = 0; 
switch model
    case 'flexible'
        isviolated = sum(x(:,end-(length(exppriorityVec)-2):end),2) >= 1; % violates if p_high + p_med >= 1
    case 'min_error'
        isviolated = isviolated | (exp(x(:,1))./exp(x(:,2))) <= 3.*(exp(x(:,end))./2); % k > 3*gamma/2
        isviolated = isviolated | (exp(x(:,end)).*3 > exp(x(:,1))); % Jbar_total > gamma*3
        isviolated = isviolated | (exp(x(:,1)) <= 3.*exp(x(:,2))); % Jbar_total > 3*tau
end

end

