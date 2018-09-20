function isviolated = check_nonbcon(model,x)
%CHECK_NONBCON checks for any violations in non-bound constraints

isviolated = 0; 
switch model
    case 'flexible'
        isviolated = sum(x(:,end-1:end),2) >= 1; % violates if p_high + p_med >= 1
    case 'min_error'
        isviolated = isviolated | (exp(x(:,1))./exp(x(:,2))) <= 3.*(exp(x(:,end))./2); % k > 3*psi/2
        isviolated = isviolated | (exp(x(:,end)).*3 > exp(x(:,1))); % Jbar_total > psi*3
        isviolated = isviolated | (exp(x(:,1)) <= 3.*exp(x(:,2))); % Jbar_total > 3*tau
end

end

