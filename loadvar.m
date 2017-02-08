function varargout = loadvar(varargin)

nvars = length(varargin);
varargout = cell(1,nvars);
for ivar = 1:nvars;
    var = varargin{ivar};
    switch var
        case 'JVec'
            nJSamp = 100;
            JVec = linspace(1e-5,20,nJSamp); % ASPEN: make sure range is reasonable for parameter range
            varargout{ivar} = JVec;
        case 'rVec'
            % radius stuff
            nRs = 200;
            rVec = linspace(0,5,nRs); % ASPEN: make sure this range is reasonable
            varargout{ivar} = rVec;
    end
end