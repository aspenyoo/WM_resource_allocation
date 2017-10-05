function countt = nonbcon(model,x)

countt = 0; 
switch model
    case 2
        countt = sum(x(:,end-1:end),2) >= 1; % violates if p_high + p_med >= 1
    case 4
        countt = countt | (exp(x(:,1))./exp(x(:,2))) <= 3.*(exp(x(:,end))./2); % k > 3*psi/2
        countt = countt | (exp(x(:,end)).*3 > exp(x(:,1))); % Jbar_total > psi*3
        countt = countt | (exp(x(:,1)) <= 3.*exp(x(:,2))); % Jbar_total > 3*tau
end
end