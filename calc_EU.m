% calculate expected utility 



Jbar_total = 10;
qVec = [0.1 0.3 0.6]; % from low to high
assert(sum(qVec) == 1)

JbarVec = qVec*Jbar_total;

