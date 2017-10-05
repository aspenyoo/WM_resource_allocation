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