function barPlotTaskFix(r2,Size)

% mybar(0.4+rand(2,3),'groupLabels',{'Group1','Group2'},'withinGroupLabels',{'Value1','Value2','Value3'},...
%     'yAxisMin=0','yAxisMax=1.5','yError',rand(2,3)*0.1,'withinGroupColors',{[0.3 0.7 0.2] [0.4 0.2 0.1] [0.1 0.0 0.0]})

early = {'V1';'V2';'V3';'V4'};
IPS = {'IPS0'; 'IPS1';'IPS2';'IPS3';'ISP4'};

cond = {'Task','Stimulus'};
% brewer = brewermap(6,'Spectral');
% condColors = {brewer(6,:) brewer(5,:)};
% condColors= {[0 0.6 0.8] [0.64 0.8 0.35]};
% condColors= {[0.2 0.51 0.73] [0.65 0.85 0.41]};

brewer = brewermap(9,'YlGnBu');
condColors = {brewer(6,:) [.65 .87 .90]};
% condColors = {brewer(6,:) brewer(2,:)};

sizedata = Size(1:4,2:3);
sizeerr = Size(1:4,4:5);

r2data = r2(1:4,2:3);
% type='grouped';

 
%  ax = gca;
%  retval.barHandles = bar(ax,data,type,'ShowBaseLine','off','BaseValue',baseValue);%, 'barWidth',0.85);
%  
%  
%   for i = 1:nBarsInGroup
%     set(retval.barHandles(i),'FaceColor',condColors{i});
%     set(retval.barHandles(i),'EdgeColor',[0.7,0.7,0.7]); %[0 0 0]'None'
%   end
%   
earlyLabel = {'V1';'V2';'V3';'hV4'};
mybar2(r2(1:4,2:3), r2(1:4,4:5), 'r^2', condColors, earlyLabel, 0.65)
mybar2(Size(1:4,2:3), Size(1:4,4:5), 'pRF size (deg)', condColors, earlyLabel, 0.65)

% mybar2(r2(1:4,2:3), r2(1:4,4:5), 'r^2', condColors, early)
% mybar2(Size(1:4,2:3), Size(1:4,4:5), 'pRF size', condColors, early)

mybar2(r2(7:11,2:3), r2(7:11,4:5), 'r^2', condColors, IPS, 0.65)
mybar2(Size(7:11,2:3), Size(7:11,4:5), 'pRF size (deg)', condColors, IPS,0.65)

% mybar2(r2(8:12,2:3), r2(8:12,4:5), 'r^2', condColors, early)
% mybar2(Size(8:12,2:3), Size(8:12,4:5), 'pRF size', condColors, early)

  
% function h=BarSpecial(data, overallWidth )
%     [r,c] = size(data);

function retval = mybar2(vars, err, yLab, condColors, xTickLab, overallWidth)
figure
baseValue = 0;
nGroups = size(vars,1);
nBarsInGroup = size(vars,2);

    retval.barHandels = zeros(nBarsInGroup,1);
    width = overallWidth / nBarsInGroup;
    offset = [-width/2 width/2];
    for i=1:nBarsInGroup
        retval.barHandles(i) = bar(vars(:,i),'FaceColor',condColors{i},'BarWidth',width-.025, 'EdgeColor',[0.6 0.6 0.6]);   
        set(retval.barHandles(i),'XData',get(retval.barHandles(i),'XData')+offset(i));
        hold on               
    end    
% end
ax = gca;
set(ax, 'XTick',1:length(xTickLab));  
set(gca,'XTickLabel', xTickLab);


for i = 1:nBarsInGroup
    barXpos(i,:) = get(retval.barHandles(i),'XData');
end
 
 
 if ~isempty(err)
  if size(err,1) ~= nGroups
    disp(sprintf('(mybar) yError does not have enough groups (%i, should be %i)',size(yError,1),nGroups));
  elseif size(err,2) ~= nBarsInGroup
    disp(sprintf('(mybar) yError does not have enough within groups (%i, should be %i)',size(yError,2),nBarsInGroup));
  else
    for iWithinGroup = 1:nBarsInGroup
      for iGroup = 1:nGroups
	% error plots up for positive values, down for negative values
        if vars(iGroup,iWithinGroup) >= baseValue
             thisYerror = err(iGroup,iWithinGroup);
        else
             thisYerror = -err(iGroup,iWithinGroup);
        end
	% plot the line
	retval.yErrorHandle(iGroup,iWithinGroup) = line([barXpos(iWithinGroup,iGroup) barXpos(iWithinGroup,iGroup)],[vars(iGroup,iWithinGroup)-thisYerror vars(iGroup,iWithinGroup)+thisYerror],'Color',[0 0 0]);
%     	retval.yErrorHandle(iGroup,iWithinGroup) = line([barXpos(iWithinGroup,iGroup) barXpos(iWithinGroup,iGroup)],[vars(iGroup,iWithinGroup) vars(iGroup,iWithinGroup)-thisYerror],'Color',[0 0 0]);

      end
    end
  end
%   % get text position
%   yTextPosMax = nan(size(vars));
%   yTextPosMax(vars>=baseValue) = vars(vars>=baseValue)+err(vars>=baseValue);
%   yTextPosMax(vars<baseValue) = vars(vars<baseValue)-err(vars<baseValue);
else
  %yTextPosMax = vars;
end

set(gca,'color','none')
box off
set(gca,'fontsize',15)
ylabel(yLab, 'fontsize', 18)
%drawPublishAxis('labelFontSize=18')
%drawPublishAxis('labelFontSize=20')

%drawPublishAxis

% %%
% 
% % get axis
% a = axis;
% % get yaxis
% if isempty(yAxisMin)
%   yAxisMin = a(3);
% end
% if isempty(yAxisMax)
%   yAxisMax = a(4);
% end
% yTickSize = (yAxisMax-yAxisMin)/15;
% xAxisMin = a(1);
% xAxisMax = a(2);
% xTickSize = (xAxisMax-xAxisMin)/15;
% 
% yaxis(yAxisMin,yAxisMax);
% 
%  
%  %%
% 
% %h = bar(r2(1:4,2:3))
% retval= mybar(data, 'groupLabels', early, 'withinGroupLabels',cond, 'dispValues',0, 'yLabelText','r^2', ...
%     'yError', err, 'withinGroupColors', condColors)
% 
% [4 153 203] [161 204 91]