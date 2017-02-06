function [p_Hit] = calc_p_Hit(r,sigma)
% calculates the probability that the saccade target lands within the disk
% 
% R: radius of disk
% SIGMA: memory noise

% check to make sure r is a orizontal hvector
r = r(:)'; 

edgee = 2;
nsamps = 100;
rangee = linspace(-edgee,edgee,nsamps);

[xx,yy] = meshgrid(rangee,rangee);
xx = xx(:); yy = yy(:);

% tic;
p_XgivenS = mvnpdf([xx yy],0,[sigma 0; 0 sigma]);
p_XgivenS = exp(log(p_XgivenS)- log(sum(p_XgivenS))); % normalize

idxs = bsxfun(@(x,y) x <= y, xx.^2+yy.^2, r.^2);
p_Hit = sum(bsxfun(@times,p_XgivenS,idxs));

% idxs = find(idxs);
% idxs = mod(idxs,nsamps^2);
% idxs(idxs==0) = nsamps^2; 
% pHit = sum(pXgivenS(idxs));
% % toc

% tic;
% rVec = r;
% nR = length(rVec);
% for iR = 1:nR
%     r = rVec(iR);
%     
%     pXgivenS = mvnpdf([xx yy],0,[sigma 0; 0 sigma]);
%     pXgivenS = exp(log(pXgivenS)- log(sum(pXgivenS))); % normalize
%     
%     idxs = xx.^2+yy.^2 <= r.^2;
%     pHit = sum(pXgivenS(idxs));
% end
% toc