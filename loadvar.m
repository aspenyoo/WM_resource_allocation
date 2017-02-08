function varargout = loadvar(var)

switch var
    case 'JVec'
        nJSamp = 100;
        JVec = linspace(1e-5,20,nJSamp); % ASPEN: make sure range is reasonable for parameter range
        varargout = {JVec};

end