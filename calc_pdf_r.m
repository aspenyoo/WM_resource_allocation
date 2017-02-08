function pdf_r = calc_pdf_r(beta, JVec)
%
% calculates the probilility density function of r for a given beta and J (or vector
% JVec)

rVec = loadvar('rVec');

% EU(r,J) for each combination of r and J in rVec and jVec
EU = bsxfun(@times, calc_p_Hit(rVec,JVec), rewardFn(rVec)');

% p(r)
pdf_r = exp(beta.*EU);
pdf_r = bsxfun(@rdivide, pdf_r, sum(pdf_r)); % normalize across rVec