function simtheta = gensimtheta(thetaMat,lb,ub,nSamps)
if nargin < 4; nSamps = 1; end

nParams = size(thetaMat,2);

for iparam = 1:nParams
    
    bandwidthh = 0.4;
    distpdf = fitdist(thetaMat(:,iparam),'Kernel','BandWidth',bandwidthh);
    x = linspace(lb(iparam),ub(iparam),100);
    pdff = pdf(distpdf,x(2:end));
    cdff = [0 cumsum(pdff)];
    cdff = cdff./cdff(end);
    
    % inverse sampling
    unifsamps = mat2cell(rand(1,nSamps),1,ones(1,nSamps));
    idx1 = cellfun(@(x) find(x >= cdff ,1,'last'),unifsamps,'UniformOutput',false);
    idx2 = cellfun(@(x) x+1,idx1,'UniformOutput',false);
    
    propdiff = cellfun(@(u,x,y) (u-cdff(y))/(cdff(x)-cdff(y)),unifsamps,idx2,idx1,'UniformOutput',false);
    
    diffx = x(2)-x(1); % difference of intervals in x
    
    xx = cellfun(@(z,y) x(z) + diffx*y,idx1,propdiff,'UniformOutput',false);
    simtheta(:,iparam) = cell2mat(xx);
    
    % plot(x,ySix,'k-','LineWidth',2)
end


%% stuff that i copied in afterward!!

% load bound and nonbound constraints and logflag
[logflag, lb, ub] = loadconstraints(model,expnumber);

filepath = ['fits/exp' num2str(expnumber) '/'];
load([filepath 'fits_model' num2str(model) '.mat'],'ML_parameters')
ML_parameters(:,logflag) = log(ML_parameters(:,logflag));

constantt = 10;
while size(simtheta,1) < nSamps
    newTheta = gensimtheta(ML_parameters,lb,ub,constantt);
    newTheta(:,logflag) = exp(newTheta(:,logflag));
    countt = nonbcon(newTheta);
    newTheta(logical(countt),:) = [];
    simtheta = [simtheta; newTheta];
    constantt = nSubj - size(simtheta,1);
end

function countt = nonbcon(x)

countt = 0;
countt = countt | (exp(x(:,1)) <= 3.*exp(x(:,2))); % Jbar_total > 3*tau

if model == 2;
    countt = countt | (sum(x(:,end-1:end),2) >= 1);
end

if model == 4;
    countt = countt | (exp(x(:,1))./exp(x(:,2))) <= 3.*(exp(x(:,end))./2); % k > 3*psi/2
    countt = countt | (exp(x(:,end)).*3 > exp(x(:,1))); % Jbar_total > psi*3
end

end