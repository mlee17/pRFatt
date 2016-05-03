%%%%%%%%%%%%%%%%%%%%
%% ScatterROI %%
%%%%%%%%%%%%%%%%%%%%
function scatterROI(roi,varargin)%thisr2, thisX, thisY, thisEccentricity, thisRfHalfWidth, r2cutoff)
r2=1;
xpos=0;
ypos=0;
ecc=0;
width=1;

%ROIs = {'V1','V2v','V2d','V3v','V3d','V3A','V3B','LO1','LO2','V4','MT','IPS0','IPS1','IPS2','IPS3','IPS4','IPS5'};
roiList = {'V1','V2v','V2d','V3v','V3d','V3A','V3B','LO1','LO2','V4','IPS0','IPS1','IPS2','IPS3','IPS4','IPS5'};
load('bothhemi_IPS5.mat')
r2cutoff = .2;

if ~(strcmp(roi,'V2') || strcmp(roi,'V3'))
        
        roiInd = find(strcmp(roiList, roi));

%mlrSmartfig('pRFPlotMultipleROIs','reuse'); clf;
 
    thisr2 = bothHemi.thisr2{roiInd}; thisX = bothHemi.thisX{roiInd}; thisY = bothHemi.thisY{roiInd}; 
    thisEccentricity= bothHemi.thisEccentricity{roiInd}; thisRfHalfWidth=bothHemi.thisRfHalfWidth{roiInd};   
else
    
     if strcmp(roi,'V2')
            ind1 = 2; ind2 = 3;
     elseif strcmp(roi,'V3')
            ind1 = 4; ind2 = 5;
     end
     
     for cond = 1:2
            thisr2{cond} = [bothHemi.thisr2{ind1}{cond}, bothHemi.thisr2{ind2}{cond}];
            thisX{cond} = [bothHemi.thisX{ind1}{cond}, bothHemi.thisX{ind2}{cond}];
            thisY{cond} = [bothHemi.thisY{ind1}{cond}, bothHemi.thisY{ind2}{cond}];
            thisEccentricity{cond} = [bothHemi.thisEccentricity{ind1}{cond}, bothHemi.thisEccentricity{ind2}{cond}];
            thisRfHalfWidth{cond} = [bothHemi.thisRfHalfWidth{ind1}{cond}, bothHemi.thisRfHalfWidth{ind2}{cond}];
     end
end
c = [1 1 1];    
cc = [0.8 0.973 0.965];

markersize = 5.3;
minX = -37; maxX=37; minY=-22; maxY=22;
figuresize = [100 500 250 250];
     
% for s = 1:2
%     minX(s) = min(d{s}.stimX(:)); %-37.4837521796327
%     maxX(s) = max(d{s}.stimX(:)); %37.4837521796327
%     minY(s) = min(d{s}.stimY(:)); %-22.2014869145184
%     maxY(s) = max(d{s}.stimY(:)); %22.2014869145184
% end
     
if r2 ==1;
% plot r^2
        figure('Position', figuresize)
        for i = 1:length(thisr2{2})
%              h1 = plot([0 0.25 0.5 0.75],nocatch(1,:),'o','MarkerSize',15);
%     set(h1(1),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',map(1,:),'LineWidth',1.5);
                plot(thisr2{2}(i),thisr2{1}(i),'o','MarkerFaceColor',(1-c*mean([thisr2{1}(i) thisr2{2}(i)])/max(mean([thisr2{1}; thisr2{2}]))).*cc, 'MarkerEdgeColor',[1 1 1],'LineWidth', 0.3,'markersize', markersize);
                hold on;
        end
%     plot(thisr2{2}, thisr2{1}, 'k.', 'markersize', 10); hold on
%     lsline;
    %x = 0:.005:1; y=x;
    x=min([thisr2{2} thisr2{1}]): .01: 0.8; y=x;
    plot(x,y,'--', 'Color', [0.6,0.6,0.6], 'LineWidth',1);
    
    x = thisr2{2};
%     x = [x(:) ones(size(x(:)))];
    y = thisr2{1};
%     beta = x\y; %((x'*x)^-1)*(x')*y'; beta = ((x'*w*x)^-1)*(x'*w)*y';
%     if ~isempty(beta)
%         plot([0 1], [0 1]*beta(1)+beta(2), 'r-');
%     end
    if length(x)>2
        mdl =fitlm(x,y);
        %plot([0 1], [0 1]*double(mdl.Coefficients(2,1))+double(mdl.Coefficients(1,1)), '-', 'Color', [0.9 0.3 0.6]);
        %title(sprintf('%s r^2\nslope:%0.2f (%s)\noffset:%0.2f (%s) <r2=%0.2f>', roi, double(mdl.Coefficients(2,1)), pvaldisp(double(mdl.Coefficients(2,4))), ...
        %double(mdl.Coefficients(1,1)), pvaldisp(double(mdl.Coefficients(1,4))), mdl.Rsquared.Ordinary));
    else
        title(sprintf('%s r^2',roi));
    end
    maxr2 = max([thisr2{2} thisr2{1}]);
    xaxis(r2cutoff, round(maxr2*10)/10);
    yaxis(r2cutoff, round(maxr2*10)/10);
    axis square;
    xlabel('Stimulus');
    ylabel('Task');
    box off;
    %set(gca,'YTick',0.2:.2:.8); set(gca,'XTick',0.2:.2:.8);    
    ylim = get(gca,'YLim'); xlim = get(gca,'XLim');
    text(xlim(2),ylim(1),sprintf('N=%s',num2str(length(thisr2{2}))),'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',12);

end
weights = mean([thisr2{1}; thisr2{2}]);

if xpos==1;
 % plot RF center (x)
figure('Position', figuresize)
    for i = 1:length(thisX{2})
%       if ~isnan(thisr2{2}(i)) && ~isnan(thisr2{1}(i))
        plot(thisX{2}(i),thisX{1}(i),'k.','Color',(1-c*mean([thisr2{1}(i) thisr2{2}(i)])/max(mean([thisr2{1}; thisr2{2}]))).*cc, 'markersize', markersize);
        hold on;
%       end
    end
%     plot(thisX{2}, thisX{1}, 'k.', 'markersize', 10); hold on
%     lsline;
    x = minX:.1:maxX; y=x;
    plot(x,y,'--', 'Color', [0.6,0.6,0.6], 'LineWidth',1);
        x = thisX{2};
        y = thisX{1};

    if length(x)>2
        %mdl =fitlm(x,y);
        %plot([minX maxX], [minX maxX]*double(mdl.Coefficients(2,1))+double(mdl.Coefficients(1,1)), ':', 'Color', [0.9 0.6 0.7]);

        mdl =fitlm(x,y,'linear','Weights',weights);
        %plot([minX maxX], [minX maxX]*double(mdl.Coefficients(2,1))+double(mdl.Coefficients(1,1)), '-', 'Color', [0.9 0.3 0.6]);
        title(sprintf('%s RF center (x)\nslope:%0.2f (%s)\noffset:%0.2f (%s) <r2=%0.2f>', roi, double(mdl.Coefficients(2,1)), pvaldisp(double(mdl.Coefficients(2,4))), ...
        double(mdl.Coefficients(1,1)), pvaldisp(double(mdl.Coefficients(1,4))), mdl.Rsquared.Ordinary));
    else
        title(sprintf('%s RF center (x)',roi));
    end
    
    xaxis(minX, maxX);
    yaxis(minX, maxX);
    axis square;
    xlabel('Stimulus');
    ylabel('Task');
    box off; set(gca, 'XTick',[-30:10:30]); set(gca,'YTick', [-30:10:30]);
%     title(sprintf('%s RF center (x)',[hemi{h} ROIs{roi}]));
end
if ypos==1 
 % plot RF center (y)
figure('Position', figuresize)
    for i = 1:length(thisY{2})
%       if ~isnan(thisr2{2}(i)) && ~isnan(thisr2{1}(i))
        plot(thisY{2}(i),thisY{1}(i),'k.','Color',(1-c*mean([thisr2{1}(i) thisr2{2}(i)])/max(mean([thisr2{1}; thisr2{2}]))).*cc, 'markersize', markersize);
        hold on;
%       end
    end
%     plot(thisY{2}, thisY{1}, 'k.', 'markersize', 10); hold on
%     lsline;
    x = minY:.1:maxY; y=x;
    plot(x,y,'--', 'Color', [0.6,0.6,0.6], 'LineWidth',1);
    
    x = thisY{2};
        y = thisY{1};

    if length(x)>2
        %mdl =fitlm(x,y);
        %plot([minY maxY], [minY maxY]*double(mdl.Coefficients(2,1))+double(mdl.Coefficients(1,1)), ':', 'Color', [0.9 0.6 0.7]);

        mdl =fitlm(x,y,'linear','Weights',weights);
        %plot([minY maxY], [minY maxY]*double(mdl.Coefficients(2,1))+double(mdl.Coefficients(1,1)), '-', 'Color', [0.9 0.3 0.6]);
        title(sprintf('%s RF center (y)\nslope:%0.2f (%s)\noffset:%0.2f (%s) <r2=%0.2f>', roi, double(mdl.Coefficients(2,1)), pvaldisp(double(mdl.Coefficients(2,4))), ...
        double(mdl.Coefficients(1,1)), pvaldisp(double(mdl.Coefficients(1,4))), mdl.Rsquared.Ordinary));
    else
        title(sprintf('%s RF center (y)',roi));
    end
    xaxis(minY, maxY);
    yaxis(minY, maxY);
    axis square;
    xlabel('Stimulus');
    ylabel('Task');
    box off; set(gca, 'XTick',[-20:10:20]); set(gca,'YTick', [-20:10:20]);
%     title(sprintf('%s RF center (y)',[hemi{h} ROIs{roi}]));
end

if width==1
 % plot rfHalfWidth
figure('Position', figuresize)
    for i = 1:length(thisRfHalfWidth{2})
%       if ~isnan(thisr2{2}(i)) && ~isnan(thisr2{1}(i))
        plot(thisRfHalfWidth{2}(i),thisRfHalfWidth{1}(i),'o','MarkerFaceColor',(1-c*mean([thisr2{1}(i) thisr2{2}(i)])/max(mean([thisr2{1}; thisr2{2}]))).*cc, 'MarkerEdgeColor',[1 1 1],'LineWidth', 0.3,'markersize', markersize);
        hold on;
%       end
    end
%     plot(thisRfHalfWidth{2}, thisRfHalfWidth{1}, 'k.', 'markersize', 10); hold on
%     lsline
    x = min([thisRfHalfWidth{2} thisRfHalfWidth{1}]):.1:max([thisRfHalfWidth{2} thisRfHalfWidth{1}]); y=x;
    plot(x,y,'--', 'Color', [0.6,0.6,0.6], 'LineWidth',1);
    
    x = thisRfHalfWidth{2};
    y = thisRfHalfWidth{1};

    if length(x)>2
         %mdl =fitlm(x,y);
       % plot([0 max([thisRfHalfWidth{1} thisRfHalfWidth{2}])], [0 max([thisRfHalfWidth{1} thisRfHalfWidth{2}])]*double(mdl.Coefficients(2,1))+double(mdl.Coefficients(1,1)), ':', 'Color', [0.9 0.6 0.7]);

        mdl =fitlm(x,y,'linear','Weights',weights);
       % plot([0 max([thisRfHalfWidth{1} thisRfHalfWidth{2}])], [0 max([thisRfHalfWidth{1} thisRfHalfWidth{2}])]*double(mdl.Coefficients(2,1))+double(mdl.Coefficients(1,1)), '-', 'Color', [0.9 0.3 0.6]);
        %title(sprintf('%s RF half width\nslope:%0.2f (%s)\noffset:%0.2f (%s) <r2=%0.2f>', roi, double(mdl.Coefficients(2,1)), pvaldisp(double(mdl.Coefficients(2,4))), ...
        %double(mdl.Coefficients(1,1)), pvaldisp(double(mdl.Coefficients(1,4))), mdl.Rsquared.Ordinary));
    else
        title(sprintf('%s RF half width', roi));
    end
    
%         xaxis(0, max([thisRfHalfWidth{1} thisRfHalfWidth{2}]));
%         yaxis(0, max([thisRfHalfWidth{1} thisRfHalfWidth{2}]));
     xaxis(min([thisRfHalfWidth{1} thisRfHalfWidth{2}]), max([thisRfHalfWidth{1} thisRfHalfWidth{2}]));
     yaxis(min([thisRfHalfWidth{1} thisRfHalfWidth{2}]), max([thisRfHalfWidth{1} thisRfHalfWidth{2}]));
       ylim = get(gca,'YLim'); xlim = get(gca,'XLim');
    text(xlim(2),0,sprintf('N=%s',num2str(length(thisr2{2}))),'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',12);
    
%      xaxis(0,20)
%      yaxis(0,20)
    axis square;
    xlabel('Stimulus');
    ylabel('Task');
    box off;
end
   
if ecc==1
 % plot eccentricity
figure('Position', figuresize)
    for i = 1:length(thisEccentricity{2})
%       if ~isnan(thisr2{2}(i)) && ~isnan(thisr2{1}(i))
        plot(thisEccentricity{2}(i),thisEccentricity{1}(i),'k.','Color',(1-c*mean([thisr2{1}(i) thisr2{2}(i)])/max(mean([thisr2{1}; thisr2{2}]))).*cc, 'markersize', markersize);
        hold on;
%       end
    end
%     plot(thisEccentricity{2}, thisEccentricity{1}, 'k.', 'markersize', 10); hold on
%     lsline
    x = min([thisEccentricity{2} thisEccentricity{1}]):1:max([thisEccentricity{2} thisEccentricity{1}]); y=x;
    plot(x,y,'--', 'Color', [0.6,0.6,0.6], 'LineWidth',1);
    x = thisEccentricity{2};
    y = thisEccentricity{1};

    if length(x)>2
       % mdl = fitlm(x,y);
       % plot([0 max([thisEccentricity{1} thisEccentricity{2}])], [0 max([thisEccentricity{1} thisEccentricity{2}])]*double(mdl.Coefficients(2,1))+double(mdl.Coefficients(1,1)), ':', 'Color', [0.9 0.6 0.7]);
        mdl =fitlm(x,y,'linear','Weights',weights);
       % plot([0 max([thisEccentricity{1} thisEccentricity{2}])], [0 max([thisEccentricity{1} thisEccentricity{2}])]*double(mdl.Coefficients(2,1))+double(mdl.Coefficients(1,1)), '-', 'Color', [0.9 0.3 0.6]);
        title(sprintf('%s Eccentricity\nslope:%0.2f (%s)\noffset:%0.2f (%s) <r2=%0.2f>', roi, double(mdl.Coefficients(2,1)), pvaldisp(double(mdl.Coefficients(2,4))), ...
        double(mdl.Coefficients(1,1)), pvaldisp(double(mdl.Coefficients(1,4))), mdl.Rsquared.Ordinary));
    else
        title(sprintf('%s Eccentricity',roi));
    end
    
    xaxis(min([thisEccentricity{1} thisEccentricity{2}]), max([thisEccentricity{1} thisEccentricity{2}]));
    yaxis(min([thisEccentricity{1} thisEccentricity{2}]), max([thisEccentricity{1} thisEccentricity{2}]));
%      xaxis(0,50)
%      yaxis(0,50)
    axis square;
    xlabel('Stimulus');
    ylabel('Task');
    box off;
end



