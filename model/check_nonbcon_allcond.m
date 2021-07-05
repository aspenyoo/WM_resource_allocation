function isviolated = check_nonbcon_allcond(model,Theta,exppriorityVec)

% hypothesis is that tms to sPCS decreases priority effect
%                        to IPS2 lowers overall Jbartotal

model = model(5:end);

% noTMS condition
x = Theta([1 3 4]);
isviolated = check_nonbcon(model,x,exppriorityVec);

% l IPS2
x = Theta([2 3 4]);
iv = check_nonbcon(model,x,exppriorityVec);
isviolated = isviolated | iv;

% l SPCS
x = Theta([1 3 5]);
iv = check_nonbcon(model,x,exppriorityVec);
isviolated = isviolated | iv;