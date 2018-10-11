function [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)

% tight_subplot creates "subplot" axes with adjustable gaps and margins
%
% [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
%
%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width 
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins 
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins 
%
%  out:  ha     array of handles of the axes objects
%                   starting from upper left corner, going row-wise as in
%                   subplot
%        pos    positions of the axes objects
%
%  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

% Pekka Kumpulainen 21.5.2012   @tut.fi
% Tampere University of Technology / Automation Science and Engineering


if nargin<3; gap = .02; end
if nargin<4 || isempty(marg_h); marg_h = .05; end
if nargin<5; marg_w = .05; end

if numel(gap)==1
    gaps{1} = gap*ones(1,Nh-1);
    gaps{2}= gap*ones(1,Nw-1); % equal gaps down and across
elseif (~iscell(gap)) && (numel(gap) == 2)
    gaps{1} = gap(1)*ones(1,Nh-1);
    gaps{2} = gap(2)*ones(1,Nw-1);
else
    gaps = gap;
end

if numel(marg_w)==1
    marg_w = [marg_w marg_w];
end
if numel(marg_h)==1
    marg_h = [marg_h marg_h];
end

axh = (1-sum(marg_h)-sum(gaps{1}))/Nh; % height of one plot
axw = (1-sum(marg_w)-sum(gaps{2}))/Nw; % width of one plot

py = 1-marg_h(2)-axh; % y-location of first x-axis

% ha = zeros(Nh*Nw,1);
ii = 0;
for ih = 1:Nh
    px = marg_w(1); % normalized x-coordinate of first y-axis
    
    for iw = 1:Nw
        ii = ii+1;
        ha(ii) = axes('Units','normalized', ...
            'Position',[px py axw axh], ...
            'XTickLabel','', ...
            'YTickLabel','');
        if iw < Nw; px = px+axw+gaps{2}(iw); end % update x-coordinate of next y-axis
    end
    if ih < Nh; py = py-axh-gaps{1}(ih); end % update y-location of next x-axis
end
if nargout > 1
    pos = get(ha,'Position');
end
ha = ha(:);
