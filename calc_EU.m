function EU = calc_EU(rVec,JVec)

EU = bsxfun(@times, calc_p_Hit(rVec,JVec), rewardFn(rVec)');