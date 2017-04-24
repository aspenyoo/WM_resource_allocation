function EU = calc_EU(rVec,JVec,alpha)

EU = bsxfun(@times, calc_p_Hit(rVec,JVec), rewardFn(rVec,alpha)');