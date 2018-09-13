
%% supplementary: how proportion allocation changes with Jbar

clear all
% close all
exppriorityVec = [0.6 0.3 0.1];

% SET UP TERNARY PLOT
colorMat = [1 0 0; 0 0 1; 0 0 0];

% axis
figure;
[h,hg,htick]=terplot;
c1 = [1 0.5 1/3];
c2 = [0 0.5 1/3];
c3 = [0 0 1/3];
hlabels=terlabel('high','medium','low');
set(h,'LineWidth',1)

% define colors
axis1 = colorMat(2,:);
axis2 = colorMat(1,:);
axis3 = colorMat(3,:);
grey = 0.7*ones(1,3);

% change the color of the grid lines
set(hg(:,1),'color',axis1)
set(hg(:,2),'color',axis2)
set(hg(:,3),'color',axis3)

% modify the label size and color
set(hlabels,'fontsize',12)
set(hlabels(1),'color',axis1)
set(hlabels(2),'color',axis2)
set(hlabels(3),'color',axis3)

% modify the tick label colors
set(htick(:,1),'color',axis1,'linewidth',3)
set(htick(:,2),'color',axis2,'linewidth',3)
set(htick(:,3),'color',axis3,'linewidth',3)

% plot proportional model
h = ternaryc(0.6,0.3,0.1);
propcolor = [0.1 0.7 0.3];
set(h,'marker','.','markerfacecolor',propcolor,'markersize',24,'markeredgecolor',propcolor)
hold on

% CALCULATE AND PLOT MIN-ERROR AND MAX-PTS MODEL PREDICTIONS
nSamps = 10;

Jbar = 2;
tau = 0.4;
alpha = 1; 
beta = 1; 
gamma = 1;

% Vec1 = linspace(0.01,3,nSamps); % tau
Vec2 = linspace(3,25,nSamps); % multiplier to use for Jbar


[pVec_MP, pVec_ME] = deal(nan(nSamps,3));
for isamp2 = 1:nSamps
    Jbar = Vec2(isamp2)*tau + 0.01
    
    MPtheta = [Jbar tau alpha beta];
    MEtheta = [MPtheta gamma];
    
    % maximizing points
    pVec_MP(isamp2,:) = calc_pVec_maxpoints(MPtheta,exppriorityVec);
    
    % minimizing error
    pVec_ME(isamp2,:) = calc_pVec_minerror(MEtheta,exppriorityVec);
end


sz = 40;
c = linspace(1,nSamps,nSamps);
% PLOT MAX POINTS MODEL
x=0.5-pVec_MP(:,1).*cos(pi/3)+pVec_MP(:,2)./2;
y=0.866-pVec_MP(:,1).*sin(pi/3)-pVec_MP(:,2).*cot(pi/6)/2;
plot(x,y,'k-')
colormap('autumn')
s = scatter(x,y,sz,c,'filled')

% PLOT MIN ERROR MODEL
x=0.5-pVec_ME(:,1).*cos(pi/3)+pVec_ME(:,2)./2;
y=0.866-pVec_ME(:,1).*sin(pi/3)-pVec_ME(:,2).*cot(pi/6)/2;
plot(x,y,'k-')
colormap('winter')
s2 = scatter(x,y,sz,c,'filled')


%% supplementary: how proportion allocation changes with probe probability
% for ME model

clear all
close all
clc

% theta
theta = [2 0.4 1 1 1];
% ------ SET UP TERNARY PLOT ----- 
colorMat = [1 0 0; 0 0 1; 0 0 0];

% axis
figure;
[h,hg,htick]=terplot;
c1 = [1 0.5 1/3];
c2 = [0 0.5 1/3];
c3 = [0 0 1/3];
hlabels=terlabel('high','medium','low');
set(h,'LineWidth',1)

% define colors
axis1 = colorMat(2,:);
axis2 = colorMat(1,:);
axis3 = colorMat(3,:);
grey = 0.7*ones(1,3);

% change the color of the grid lines
set(hg(:,1),'color',axis1)
set(hg(:,2),'color',axis2)
set(hg(:,3),'color',axis3)

% modify the label size and color
set(hlabels,'fontsize',12)
set(hlabels(1),'color',axis1)
set(hlabels(2),'color',axis2)
set(hlabels(3),'color',axis3)

% modify the tick label colors
set(htick(:,1),'color',axis1,'linewidth',3)
set(htick(:,2),'color',axis2,'linewidth',3)
set(htick(:,3),'color',axis3,'linewidth',3)

% get probe probabilities that will be investigated
x = 0:0.1:1;
[X,Y] = meshgrid(x,x);
X = X(:);
Y = Y(:);
idx = (X+Y) > 1;
X(idx) = [];
Y(idx) = [];
Z = 1-X-Y;

hold on
% get x and y of actual values
x=0.5-X.*cos(pi/3)+Y./2;
y=0.866-X.*sin(pi/3)-Y.*cot(pi/6)/2;

nSamps = length(X);
pVec_ME = nan(nSamps,3);
for isamp = 1:nSamps;
    exppriorityVec = [X(isamp) Y(isamp) Z(isamp)]
    
    pVec_ME(isamp,:) = calc_pVec_minerror(theta,exppriorityVec);
end

u = 0.5-pVec_ME(:,1).*cos(pi/3)+pVec_ME(:,2)./2;
v = 0.866-pVec_ME(:,1).*sin(pi/3)-pVec_ME(:,2).*cot(pi/6)/2;

u = u-x;
v = v-y;
quiver(x,y,u,v,'k')


%% supplementary: how proportion allocation changes with probe probability
% for MP model

clear all
close all
clc

% theta
theta = [2 0.4 1 1];
% ------ SET UP TERNARY PLOT ----- 
colorMat = [1 0 0; 0 0 1; 0 0 0];

% axis
figure;
[h,hg,htick]=terplot;
c1 = [1 0.5 1/3];
c2 = [0 0.5 1/3];
c3 = [0 0 1/3];
hlabels=terlabel('high','medium','low');
set(h,'LineWidth',1)

% define colors
axis1 = colorMat(2,:);
axis2 = colorMat(1,:);
axis3 = colorMat(3,:);
grey = 0.7*ones(1,3);

% change the color of the grid lines
set(hg(:,1),'color',axis1)
set(hg(:,2),'color',axis2)
set(hg(:,3),'color',axis3)

% modify the label size and color
set(hlabels,'fontsize',12)
set(hlabels(1),'color',axis1)
set(hlabels(2),'color',axis2)
set(hlabels(3),'color',axis3)

% modify the tick label colors
set(htick(:,1),'color',axis1,'linewidth',3)
set(htick(:,2),'color',axis2,'linewidth',3)
set(htick(:,3),'color',axis3,'linewidth',3)

% get probe probabilities that will be investigated
x = 0:0.1:1;
[X,Y] = meshgrid(x,x);
X = X(:);
Y = Y(:);
idx = (X+Y) > 1;
X(idx) = [];
Y(idx) = [];
Z = 1-X-Y;


hold on
% get x and y of actual values
x=0.5-X.*cos(pi/3)+Y./2;
y=0.866-X.*sin(pi/3)-Y.*cot(pi/6)/2;

nSamps = length(X);
pVec_MP = nan(nSamps,3);
for isamp = 1:nSamps;
    exppriorityVec = [X(isamp) Y(isamp) Z(isamp)]
    
    pVec_MP(isamp,:) = calc_pVec_maxpoints(theta,exppriorityVec);
end

u = 0.5-pVec_MP(:,1).*cos(pi/3)+pVec_MP(:,2)./2;
v = 0.866-pVec_MP(:,1).*sin(pi/3)-pVec_MP(:,2).*cot(pi/6)/2;

u = u-x;
v = v-y;
quiver(x,y,u,v,'k')