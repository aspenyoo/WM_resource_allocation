function [p_Hit] = calc_p_Hit(r,J)
% calculates the probability that the saccade target lands within the disk
% 
% R: radius of disk. can be a scalar or vector
% J: memory precision. can be a scalar or vector
% 
% ============ OUTPUT VARIABLES ===========
% P_HIT: an nR x nJ vector of p(Hit) for each combination of r and J. 

r = r(:)'; % check to make sure r is a horizontal vector
J = J(:); % J is vertical


edgee = 2;
nsamps = 100;
rangee = linspace(-edgee,edgee,nsamps);

[xx,yy] = meshgrid(rangee,rangee);
xx = xx(:); yy = yy(:);

% get covariance matrix
nJs = length(J);
Sigma = zeros(1,2,nJs*nsamps^2);
Sigma(1,:,:) = sort(repmat(sqrt(1./J),nsamps^2,2),'descend')'; % in descending order to keep J ascending

% p(XgivenS)
p_XgivenS = mvnpdf(repmat([xx yy],nJs,1),0,Sigma);
p_XgivenS = reshape(p_XgivenS,nsamps^2,1,nJs);
p_XgivenS = exp(bsxfun(@minus,log(p_XgivenS),log(sum(p_XgivenS)))); % normalize

idxs = bsxfun(@(x,y) x <= y, xx.^2+yy.^2, r.^2);
p_Hit = squeeze(sum(bsxfun(@times,p_XgivenS,idxs)));