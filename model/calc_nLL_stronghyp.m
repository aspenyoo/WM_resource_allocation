function nLL = calc_nLL_stronghyp(Theta,data,exppriorityVec,fixparams,condhypVec)

load('fittingsettings.mat','condVec','nConds')
% condVec = {'noTMS','IPS2','sPCS'};

nLL = 0;

for icond = 1:nConds
    cond = condVec{icond};
    condhyp = condhypVec(icond);
    
    switch condhyp
        case 1
            theta = Theta([1 2 3]);
        case 2
            theta = Theta([1 2 5]); % diff pVec
        case 3
            theta = Theta([4 2 3]); % diff Jbar
    end
    
    nLL = nLL + calc_nLL('flexible',theta,data.(cond),exppriorityVec,fixparams);
end