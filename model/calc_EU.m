function EU = calc_EU(rVec,JVec,alpha)
% size: rVec = 1:nR

% calc_p_Hit gives size rVec x JVec
% rewardFn gives 1 x rVec

EU = bsxfun(@times, calc_p_Hit(rVec,JVec), rewardFn(rVec,alpha)');

% pHitMat = calc_p_Hit(rVec,JVec);
% EU = nan(length(rVec), length(JVec));
% for ir = 1:length(rVec);
%     EU(ir,:) = pHitMat(ir,:).*rewardFn(rVec(ir),alpha);
% end