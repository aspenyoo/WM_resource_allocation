function [AIC, BIC, AICc ] = modcomp(nLL,K,n)
% [AIC, BIC, AICc] = modcomp(nLL,K,n) computes AIC, BIC, and/or 
% AICc scores given nLL, K, AND n. 
% 
% nLL: negative log likelihood (of ML). size(nLL) = [nSubj, nModels]
% K: number of parameters. array of length nModels.
% n: number of trials. array with length nSubj.


n = n(:); % vertical
K = K(:)'; % horizontal

% checking that everything is arranged properly
assert(size(nLL,1) == length(n));
assert(size(nLL,2) == length(K));


AIC = 2*bsxfun(@plus,nLL,K);

BIC = 2*nLL + bsxfun(@times,log(n),K);

if nargout > 2
AICc = AIC + 2.*bsxfun(@rdivide,K.*(K+1),(bsxfun(@minus,n,K)-1));
end