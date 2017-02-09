

subjnum = 4;

load('cleandata.mat')
subjdata = data{subjnum};

[A,b,Aeq,beq,nonlcon] = deal([]);
options = optimset('Display','iter');
lb = [1e-5 1e-5 1e-5]; % Jbar_total, tau, beta
ub = [50 50 50]; % ASPEN: refine
logflag = logical([1 1 0]);
lb(logflag) = log(lb(logflag));
ub(logflag) = log(ub(logflag));
startingval = lb + rand(1,3).*(ub - lb);
[x,fval] = fmincon(@(x) calc_nLL(subjdata,x),startingval,A,b,Aeq,beq,lb,ub,nonlcon,options);