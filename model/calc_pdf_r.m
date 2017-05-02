function pdf_r = calc_pdf_r(beta, JVec, alpha)
%
% calculates the probilility density function of choosing r for a given beta and J (or vector
% JVec)
% 
% 
% -----------------------
%      Aspen H. Yoo
%   aspen.yoo@nyu.edu

rVec = loadvar('rVec');
deltar = diff(rVec(1:2)); % the size of the bin. this is so that the pdf is not a function ngrids

% EU(r,J) for each combination of r and J in rVec and jVec
EU = bsxfun(@times, calc_p_Hit(rVec,JVec), rewardFn(rVec,alpha)'); % SIZE: (rVec x JVec) = (rVec x JVec) x (rVec x 1)

% p(r)
pdf_r = exp(beta.*EU); % size: (rVec x JVec)
pdf_r = bsxfun(@rdivide, pdf_r, sum(pdf_r))./deltar; % normalize across rVec. size: (rVec x JVec) = (rVec x JVec) x (1 x JVec)