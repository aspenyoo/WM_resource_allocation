% bias analysis

% NOTES TO SELF
% \theta = 0 is on the positive x axis. goes from -pi to pi.
% \rho is around 10 deg.

saccadeType = 'final'; % 'primary' or 'finpal'
coordinateType = 'rt'; % xy (cartesian) or rt (rho and theta; polar)
biasType = 'axis'; % quadrant (directional bias within each quadrant) or
% axis (collapse across quadrants and see if there are biases toward/away
% from cardinal axes)
subjid = '1';
filepath = 'DATA';
isplot = [1 1];
issave = [0 0];

load([filepath '/subj' subjid '_ii_stats_targetvalue.mat']); % loading appropriate data

colorMat = aspencolors(4,'passport');

% getting primary, final, and corrective saccade in proper coordinate space
primary = []; final = []; corrective = []; condition = []; inconsistentsessions = [];
nSessions = length(ii_stats);
for isession = 1:nSessions;
    condition = [condition; ii_stats(isession).target_value];
    
    if sum(size(ii_stats(isession).target_value) == size(ii_stats(isession).primary_x)) == 2
        switch coordinateType
            case 'xy' % x-y coordinates
                primary = [primary; ii_stats(isession).primary_x ii_stats(isession).primary_y];
                final = [final; ii_stats(isession).final_x ii_stats(isession).final_y];
                corrective = [corrective; ii_stats(isession).corrective_x ii_stats(isession).corrective_y;];
                
            case 'rt' % rho-theta coordinates
                primary = [primary; ii_stats(isession).primary_rho ii_stats(isession).primary_theta];
                final = [final; ii_stats(isession).final_theta ii_stats(isession).final_rho]; % rho and theta are flipped in data
                corrective = [corrective; ii_stats(isession).corrective_theta ii_stats(isession).corrective_rho]; % rho and theta are flipped in data
        end
    else % document which sessions don't have priority info
        inconsistentsessions = [inconsistentsessions isession];
    end
end
clear ii_stats % clear large struct for speed

% get quadrant index
switch coordinateType
    case 'rt'
        quadVec = [-pi -pi/2 0 pi/2 pi]; % quads 3 4 1 2
        quadIdx = [3 4 1 2];
    case 'xy'
        quadVec = [1 1; -1 1; -1 -1; 1 -1]; % quads 1 2 3 4
        quadIdx = 1:4;
end
quadrant = nan(size(condition));
for iquad = 1:4;
    
    % get the stuff in the quadrant
    switch coordinateType
        case 'rt'
            idx = (corrective(:,2) >= quadVec(iquad)) & (corrective(:,2) < quadVec(iquad+1));
        case 'xy'
            idx = (sign(corrective(:,1)) == quadVec(iquad,1)) &  (sign(corrective(:,2)) == quadVec(iquad,2));
    end
    quadrant(idx) = quadIdx(iquad);
end

switch biasType
    case 'axis'
        switch coordinateType
            case 'xy'
                corrective = abs(corrective);
                primary(quadrant == 2,1) = -primary(quadrant == 2,1);
                primary(quadrant == 3,:) = -primary(quadrant == 3,:);
                primary(quadrant == 4,2) = -primary(quadrant == 4,2);
                final(quadrant == 2,1) = -final(quadrant == 2,1);
                final(quadrant == 3,:) = -final(quadrant == 3,:);
                final(quadrant == 4,2) = -final(quadrant == 4,2);
            case 'rt'
                primary(quadrant == 2,2) = pi - primary(quadrant == 2,2);
                primary(quadrant == 3,2) = primary(quadrant == 3,2) + pi;
                primary(quadrant == 4,2) = -primary(quadrant == 4,2);
                final(quadrant == 2,2) = pi - final(quadrant == 2,2);
                final(quadrant == 3,2) = final(quadrant == 3,2) + pi;
                final(quadrant == 4,2) = -final(quadrant == 4,2);
                corrective(quadrant == 2,2) = pi - corrective(quadrant == 2,2);
                corrective(quadrant == 3,2) = corrective(quadrant == 3,2) + pi;
                corrective(quadrant == 4,2) = -corrective(quadrant == 4,2);
        end
end

% error_primary = corrective - primary;
% error_final = corrective - final;

% % checknig to see if any changes quadrants from primary/final to corrective
% primaryaxis = primary./abs(primary);
% correctiveaxis = corrective./abs(corrective);
% idx = sum(primaryaxis - correctiveaxis,2);

% % plotting primary corrective
% plot([primary(:,1) corrective(:,1)]',[primary(:,2) corrective(:,2)]','k');
% hold on
% plot(corrective(:,1)',corrective(:,2)','k*'); % plotting final saccades
% plot([-15 15; 0 0]',[0 0; -15 15]','k');
% axis([-15 15 -15 15])
% defaultplot;
% xlabel('x-axis')
% ylabel('y-axis')

% cells of the data by condition
conditionVec = unique(condition);
nCond = length(conditionVec);
primaryCell = cell(1,nCond); finalCell = primaryCell;
idxCell = primaryCell; correctiveCell = primaryCell; quadCell = idxCell;

for icond = 1:nCond;
    primaryCell{icond} = primary(condition == conditionVec(icond),:);
    finalCell{icond} = final(condition == conditionVec(icond),:);
    correctiveCell{icond} = corrective(condition == conditionVec(icond),:);
    idxCell{icond}= idx(condition == conditionVec(icond));
    quadCell{icond} = quadrant(condition == conditionVec(icond));
end

if (isplot(1))
    % ========= PLOT "RAW" DATA (acc to biasType and coordinateType) =========
    figure;
    for icond = 1:nCond;
        switch saccadeType
            case 'primary'
                v1 = primaryCell{icond};
            case 'final'
                v1 = primaryCell{icond};
        end
        v2 = correctiveCell{icond};
        
        % plot primary/final to corrective saccade for each priority
        subplot(2,2,icond)
        hold on
        plot([v1(:,1) v2(:,1)]',[v1(:,2) v2(:,2)]','k'); % primary/final to corrective saccade
        plot(v2(:,1)',v2(:,2)','k*'); % plotting corrective saccades
        defaultplot;
        title(sprintf('priority = %0.2f',conditionVec(icond)))
        switch coordinateType % axes and stuff like that
            case 'xy'
                plot([-15 15; 0 0]',[0 0; -15 15]','Color',0.8*ones(1,3)); % axes
                axis equal
                if strcmp(biasType,'quad');
                    axis([-15 15 -15 15])
                else
                    plot([0 15/sqrt(2)],[0 15/sqrt(2)],'Color',0.8*ones(1,3)); % plotting diagonal
                    plot([0 tand(22.5)*15; 0 15*cosd(22.5)]',[0 15*cosd(22.5); 0 tand(22.5)*15]','--','Color',0.8*ones(1,3)); % plotting off axis things
                    axis([-1 15 -1 15])
                end
                xlabel('x-axis')
                ylabel('y-axis')
            case 'rt'
                xlabel('\rho')
                ylabel('\theta')
                if strcmp(biasType,'quad')
                    axis([5 15 -pi pi])
                    plot(repmat([5; 15],1,4),repmat([-pi/2 0 pi/2 pi],2,1),'Color',aspencolors(1,'blue')); % coordinate axes in polar space
                else
                    axis([5 15 -0.5 2])
                    plot(repmat([5; 15],1,5),repmat([0 pi/8 pi/4 3*pi/8 pi/2],2,1),'Color',aspencolors(1,'blue')); % pizza slices in polar space
                end
        end
        
        switch biasType
            case 'quad'
                % look at mean responses per quadrant
                meanv1 = nan(4,2); meanv2 = nan(4,2);
                for iquad = 1:4;
                    meanv1(iquad,:) = mean(v1(quadCell{icond}==iquad,:));
                    meanv2(iquad,:) = mean(v2(quadCell{icond}==iquad,:));
                end
                plot([meanv1(:,1) meanv2(:,1)]',[meanv1(:,2) meanv2(:,2)]','Color',...
                    aspencolors('coral'),'LineWidth',3)
                plot(meanv2(:,1),meanv2(:,2),'.','Color',aspencolors('coral'),'MarkerSize',18);
            case 'axis'
                meanv1 = nan(4,2); meanv2 = nan(4,2);
                switch coordinateType
                    case 'xy'
                        sliceratio = [0 tand(22.5)/cosd(22.5) 1 cosd(22.5)/tand(22.5) Inf];
                        for islice = 1:4;
                            meanv2(islice,:) = mean(v2((v2(:,2)./v2(:,1) > sliceratio(islice)) & (v2(:,2)./v2(:,1) < sliceratio(islice+1)),:));
                            meanv1(islice,:) = mean(v1((v1(:,2)./v1(:,1) > sliceratio(islice)) & (v1(:,2)./v1(:,1) < sliceratio(islice+1)),:));
                        end
                    case 'rt'
                        sliceratio = [0 pi/8 pi/4 3*pi/8 pi/2];
                        for islice = 1:4;
                            meanv2(islice,:) = mean(v2((v2(:,2) > sliceratio(islice)) & (v2(:,2) < sliceratio(islice+1)),:));
                            meanv1(islice,:) = mean(v1((v1(:,2) > sliceratio(islice)) & (v1(:,2) < sliceratio(islice+1)),:));
                        end
                end
%                 plot([meanv1(:,1) meanv2(:,1)]',[meanv1(:,2) meanv2(:,2)]','Color',...
%                     aspencolors('coral'),'LineWidth',3)
%                 plot(meanv2(:,1),meanv2(:,2),'.','Color',aspencolors('coral'),'MarkerSize',18);
                plot([0 ; 0 ]',[0 15*cosd(22.5); 0 tand(22.5)*15]','--','Color',0.8*ones(1,3))
        end
    end
end
% fig = gcf;
% saveas(fig,['imgs/subj' subjid '_' coordinateType '_' saccadeType '_quadBias.png'])
% saveas(fig,['imgs/subj' subjid '_' coordinateType '_' saccadeType '_quadBias.fig'])
%
quadIdx = [2 1 3 4];
% ====== PLOTTING ERROR DISTRIBUTIONS =======

switch saccadeType
    case 'primary'
        v1 = primary;
    case 'final'
        v1 = primary;
end
v2 = corrective;
blaherror = v2-v1;

bias = cell(1,nCond);
if (isplot(2))
    figure;
    for icond = 1:nCond;
        
        switch saccadeType
            case 'primary'
                v1 = primaryCell{icond};
            case 'final'
                v1 = primaryCell{icond};
        end
        v2 = correctiveCell{icond};
        
        
        switch biasType
            case 'quad'
                figure;
                for iquad = 1:4;
                    subplot(2,2,quadIdx(iquad)); hold on
                    vv1 = v1(quadCell{icond}==iquad,:);
                    vv2 = v2(quadCell{icond}==iquad,:);
                    error_sacc = vv2-vv1;
                    covMat = cov(error_sacc);
                    circle = bsxfun(@plus,mean(error_sacc),[cos(linspace(0,2*pi,50))' sin(linspace(0,2*pi,50))']*covMat); % calculate covariance circle
                    plot(error_sacc(:,1),error_sacc(:,2),'.k','MarkerSize',12); % error data
                    plot(mean(error_sacc(:,1)), mean(error_sacc(:,2)),'.','Color',aspencolors('coral'),'MarkerSize',12); % mean error
                    plot(circle(:,1),circle(:,2),'Color',aspencolors('coral')); % plot covariance circle
                    switch coordinateType
                        case 'xy'
                            maxerr = max(abs(error_sacc(:)));
                            plot([-maxerr maxerr; 0 0]',[0 0; -maxerr maxerr]','Color',0.8*ones(1,3)); % coordinate axes
                            axis equal
                            axis([-maxerr maxerr -maxerr maxerr])
                            xlabel('x-axis')
                            ylabel('y-axis')
                        case 'rt'
                            maxerr = max(abs(error_sacc));
                            xlabel('\rho')
                            ylabel('\theta')
                            axis([-maxerr(1) maxerr(1) -maxerr(2) maxerr(2)])
                            plot(repmat([-maxerr(1) maxerr(1)],1,1),repmat([maxerr(2) maxerr(2)],1,1),'Color',aspencolors(1,'blue')); % coordinate axes in polar space
                    end
                    defaultplot
                    title(sprintf('quadrant %d',iquad))
                    
                    biglabelplot(sprintf('priority = %0.2f',conditionVec(icond)))
                end
            case 'axis' % assuming no effect of quadrant. only effect of how close to cardinal axis
                switch coordinateType
                    case 'xy'
                        subplot(2,2,icond)
                        sliceratio = [0 tand(22.5)/cosd(22.5) 1 cosd(22.5)/tand(22.5) Inf];
                        hold on; 
                        bias{icond} = nan(4,2);
                        for islice = 1:4;
                            idx = (v2(:,2)./v2(:,1) > sliceratio(islice)) & (v2(:,2)./v2(:,1) < sliceratio(islice+1));
                            vv2 = v2(idx,:);
                            vv1 = v1(idx,:);
                            
                            % plotting errors 
                            error_sacc = bsxfun(@plus,vv2 - vv1, mean(vv2));
                            plot(error_sacc(:,1),error_sacc(:,2),'.','Color',colorMat(islice,:),'MarkerSize',12); % plotting error (centered around mean(vv1)for data visualization)
                            bias{icond}(islice,:) = mean(error_sacc) - mean(vv2);
                            plot(mean(vv2(:,1)),mean(vv2(:,2)),'*','Color',colorMat(islice,:),'MarkerSize',12); % mean(vv1)
                            
                            % plotting error distribution stuff
                            covMat = cov(error_sacc);
                            circle = bsxfun(@plus,mean(error_sacc),[cos(linspace(0,2*pi,50))' sin(linspace(0,2*pi,50))']*covMat); % calculate covariance circle
                            plot(mean(error_sacc(:,1)), mean(error_sacc(:,2)),'.','Color',aspencolors('coral'),'MarkerSize',12); % mean error
                            plot(circle(:,1),circle(:,2),'Color',aspencolors('coral')); % plot covariance circle
                            
                            % plotting relavent axes/lines
                            maxerr = max(abs(error_sacc(:)));
                            plot([-maxerr maxerr; 0 0]',[0 0; -maxerr maxerr]','Color',0.8*ones(1,3)); % coordinate axes
                            plot([0 15/sqrt(2)],[0 15/sqrt(2)],'Color',0.8*ones(1,3)); % plotting diagonal
                            plot([0 tand(22.5)*15; 0 15*cosd(22.5)]',[0 15*cosd(22.5); 0 tand(22.5)*15]','--','Color',0.8*ones(1,3)); % plotting off axis things
                            
                            
                            % plotting stuff
                            axis equal
                            axis([-15 20 -15 20])
                            %                             axis([min(error_sacc(:)) maxerr min(error_sacc(:)) maxerr])
                            %                             axis([min(v2(:) - v1(:)) max(v2(:) - v1(:)) min(v2(:) - v1(:)) max(v2(:) - v1(:))])
                            %                             axis([min(blaherror(:,1)) max(blaherror(:,1)) min(blaherror(:,2)) max(blaherror(:,2))]);
                            
                            title(sprintf('priority = %0.2f',conditionVec(icond)))
                        end
                        defaultplot
                    case 'rt'
                        subplot(2,2,icond)
                        sliceratio = [0 tand(22.5)/cosd(22.5) 1 cosd(22.5)/tand(22.5) Inf];
                        hold on; 
                        bias{icond} = nan(4,2);

                        for islice = 1:4;
                            idx = (v2(:,2) > sliceratio(islice)) & (v2(:,2) < sliceratio(islice+1));
                            idx
                            vv2 = v2(idx,:);
                            vv1 = v1(idx,:);
                            
                            % plotting errors 
                            error_sacc = bsxfun(@plus,vv2 - vv1, mean(vv2));
                            plot(error_sacc(:,1),error_sacc(:,2),'.','Color',colorMat(islice,:),'MarkerSize',12); % plotting error (centered around mean(vv1)for data visualization)
                            bias{icond}(islice,:) = mean(error_sacc) - mean(vv2);
                            plot(mean(vv2(:,1)),mean(vv2(:,2)),'*','Color',colorMat(islice,:),'MarkerSize',12); % mean(vv1)
                            
                            % plotting error distribution stuff
                            covMat = cov(error_sacc);
                            circle = bsxfun(@plus,mean(error_sacc),[cos(linspace(0,2*pi,50))' sin(linspace(0,2*pi,50))']*covMat); % calculate covariance circle
                            plot(mean(error_sacc(:,1)), mean(error_sacc(:,2)),'.','Color',aspencolors('coral'),'MarkerSize',12); % mean error
                            plot(circle(:,1),circle(:,2),'Color',aspencolors('coral')); % plot covariance circle
                            
                            
                            
                            % axis stuff
                            axis([4 16 -6 2])
                            title(sprintf('priority = %0.2f',conditionVec(icond)))
                        end
                        defaultplot
                        
                end
                
                
        end
    end
    
end


switch coordinateType % axes and stuff like that
    case 'xy'
        
    case 'rt'
        
end
