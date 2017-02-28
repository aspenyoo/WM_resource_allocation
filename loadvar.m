function varargout = loadvar(varargin)
nvars = length(varargin);
varargout = cell(1,nvars);
for ivar = 1:nvars;
    var = varargin{ivar};
    if iscell(var); Jbar = var{2}; tau = var{3}; var = var{1}; end
    switch var
        case 'JVec'
            nJSamp = 100;
            xmin = Jbar;
            xmax = Jbar;
%             if exist('Jbar')

                % lower bound
                while gampdf(xmin,Jbar/tau,tau) > 1e-4;
                    if xmin <= 1; 
                        xmin = 1e-10;
                        break
                    else
                        xmin = xmin - 1;
                    end
                end
                % upper bound
                while gampdf(xmax,Jbar/tau,tau) > 1e-4;
                    xmax = xmax + 0.1;
                end
%             end
            JVec = linspace(xmin,xmax,nJSamp);
%             JVec = linspace(1e-5,20,nJSamp); % ASPEN: make sure range is reasonable for parameter range
            varargout{ivar} = JVec;
        case 'rVec'
            % radius stuff
            nRs = 500;
            rVec = linspace(0,10,nRs);
            varargout{ivar} = rVec;
    end
end