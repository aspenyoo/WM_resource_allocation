function [p_Hit] = calc_p_Hit(r,J)
% calculates the probability that the saccade target lands within the
% circle wager. 
% 
% ====== INPUT VARIABLES =====
% R: radius of disk. can be a scalar or vector
% J: memory precision. can be a scalar or vector
% 
% ===== OUTPUT VARIABLES =====
% P_HIT: an nR x nJ vector of p(Hit) for each combination of r and J. 

r = r(:); % r is vertical
J = J(:)'; % J is horizontal

p_Hit = bsxfun(@(x,prec) 1-exp(-(prec.*(x.^2))./2),r,J);
% p_Hit = nan(length(r),length(J));
% for ir = 1:length(r);
%     p_Hit(ir,:) = 1-exp(-(r(ir).^2)./(2./J));
% end