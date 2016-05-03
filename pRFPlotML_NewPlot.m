% pRFPlot.m
%
%        $Id:$ 
%      usage: pRFPlot(v,overlayNum,scanNum,x,y,z,,roi)
%         by: justin gardner
%       date: 11/22/11
%    purpose: plot function for displaying results of pRF analysis
%
function pRFPlotML(v,overlayNum,defaultscanNum,x,y,z,roi)

% check arguments
if ~any(nargin == [7])
  help pRFPlot
  return
end

% see if the shift key is down
%shiftDown = any(strcmp(get(viewGet(v,'figureNumber'),'CurrentModifier'),'shift'));
shiftDown = any(strcmp(get(viewGet(v,'figureNumber'),'SelectionType'),'extend'));

% check if pRF has been run
a = viewGet(v,'Analysis');
if ~isfield(a,'type') || ~strcmp(a.type,'pRFAnal')
  disp(sprintf('(pRFPlot) pRF analysis has not been run on this scan'));
  return
end

% get the d
for scanNum = 1:2
d{scanNum} = viewGet(v,'d',scanNum);
if isempty(d{scanNum}),disp(sprintf('(pRFPlot) Could not find d structure for this scan %d', scanNum));return,end

% get the parametrs of the pRF fit
r2{scanNum} = viewGet(v,'overlayData',scanNum,viewGet(v,'overlayNum','r2'));
if isempty(r2{scanNum})
  disp(sprintf('(pRFPlot) pRF analysis has not been run on this scan %d', scanNum));
  return
end
thisR2(scanNum) = r2{scanNum}(x,y,z);
polarAngle{scanNum} = viewGet(v,'overlayData',scanNum,viewGet(v,'overlayNum','polarAngle'));
thisPolarAngle(scanNum) = polarAngle{scanNum}(x,y,z);
eccentricity{scanNum} = viewGet(v,'overlayData',scanNum,viewGet(v,'overlayNum','eccentricity'));
thisEccentricity(scanNum) = eccentricity{scanNum}(x,y,z);
rfHalfWidth{scanNum} = viewGet(v,'overlayData',scanNum,viewGet(v,'overlayNum','rfHalfWidth'));
thisRfHalfWidth(scanNum) = rfHalfWidth{scanNum}(x,y,z);

% get the params that have been run
scanDims{scanNum} = viewGet(v,'scanDims',scanNum);
whichVoxel(scanNum) = find(d{scanNum}.linearCoords == sub2ind(scanDims{scanNum},x,y,z));
r(scanNum) = d{scanNum}.r(whichVoxel(scanNum),:);


% if no voxel has been found in precomputed analysis then do fit (or if shift is down)
if isempty(whichVoxel(scanNum)) || shiftDown
  % check if shift is being held down, in which case we reget parameters
  if shiftDown
    fit{scanNum} = pRFFit(v,overlayNum,scanNum,x,y,z,roi);
  else
    fit{scanNum} = pRFFit(v,overlayNum,scanNum,x,y,z,roi,'fitTypeParams',a.params.pRFFit);
  end
  if isempty(fit{scanNum}),return,end
  % set the overlays
  r2{scanNum}(x,y,z) = fit{scanNum}.r2;
  polarAngle{scanNum}(x,y,z) = fit{scanNum}.polarAngle;
  eccentricity{scanNum}(x,y,z) = fit{scanNum}.eccentricity;
  rfHalfWidth{scanNum}(x,y,z) = fit{scanNum}.std;
  % reset the overlays
  v = viewSet(v,'overlayDataReplace',r2{scanNum},'r2');
  v = viewSet(v,'overlayDataReplace',polarAngle{scanNum},'polarAngle');
  v = viewSet(v,'overlayDataReplace',eccentricity{scanNum},'eccentricity');
  v = viewSet(v,'overlayDataReplace',rfHalfWidth{scanNum},'rfHalfWidth');
  % now refresh the display
  refreshMLRDisplay(viewGet(v,'viewNum'));
  return
end


params{scanNum} = d{scanNum}.params(:,whichVoxel(scanNum));
if isfield(d{scanNum},'paramsInfo')
  paramsInfo{scanNum} = d{scanNum}.paramsInfo;
else
  paramsInfo{scanNum} = [];
end


% get params
m{scanNum} = pRFFit(v,scanNum,x,y,z,'stim',d{scanNum}.stim,'getModelResponse=1','params',params{scanNum},'concatInfo',d{scanNum}.concatInfo,'fitTypeParams',a.params.pRFFit,'paramsInfo',paramsInfo{scanNum});

end

% roi
if ~shiftDown
%   pRFPlotROI(v,roi,d,a,r2,eccentricity,polarAngle,rfHalfWidth);
  pRFCompareROI(v,roi,d,a,r2,eccentricity,polarAngle,rfHalfWidth);
%   pRFCompareROI2(v,roi,d,a,r2,eccentricity,polarAngle,rfHalfWidth);
end


% and plot, set a global so that we can use the mouse to display
% different time points
global gpRFPlot;
gpRFPlot.fignum = selectGraphWin;

% clear callbacks
set(gpRFPlot.fignum,'WindowButtonMotionFcn','');
gpRFPlot.d=[]; gpRFPlot.rfModel=[];
for scanNum = 1:2
% keep the stim
gpRFPlot.d{scanNum} = d{scanNum};
gpRFPlot.rfModel{scanNum} = m{scanNum}.rfModel;
end

% keep the axis that has the time series
gpRFPlot.a = subplot(6,5,[1:4 6:9 11:14 16:19]);
% plot the rectangle that shows the current stimuli
% FIX: Start time
gpRFPlot.t = 50;
gpRFPlot.hRect = rectangle('Position',[gpRFPlot.t-4 min([m{1}.tSeries;m{2}.tSeries]) 4 max([m{1}.tSeries;m{2}.tSeries])-min([m{1}.tSeries;m{2}.tSeries])],'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
hold on
% plot time series
brewer = brewermap(9,'YlGnBu');
condColors = {brewer(6,:) [.65 .87 .90]};
plot(m{1}.tSeries,'-', 'Color', [0.3,0.7,0.8], 'Linewidth',1);%[0.7,0.5,0.6]);
plot(m{2}.tSeries,'-', 'Color', [0.5,0.5,0.5], 'Linewidth',1);%[0.5,0.7,0.8]); 
axis tight
% plot model
plot(m{1}.modelResponse,'-', 'Color', brewer(6,:),'Linewidth', 1.5);
plot(m{2}.modelResponse,'-', 'Color',[.65 .87 .90], 'Linewidth',1.5);
for scanNum = 1:2
if d{scanNum}.concatInfo.n > 1
  vline(d{scanNum}.concatInfo.runTransition(2:end,1));
end
% convert coordinates back to x,y for display
[thisx(scanNum) thisy(scanNum)] = pol2cart(thisPolarAngle(scanNum),thisEccentricity(scanNum));
end
xlabel('Time (volumes)');
ylabel('BOLD (%)');
title(sprintf('BarsTask: [%i %i %i] r^2=%0.2f polarAngle=%0.2f eccentricity=%0.2f rfHalfWidth=%0.2f %s [x=%0.2f y=%0.2f]\n%s\n BarsTaskFixation: [%i %i %i] r^2=%0.2f polarAngle=%0.2f eccentricity=%0.2f rfHalfWidth=%0.2f %s [x=%0.2f y=%0.2f]\n%s',...
    x,y,z,thisR2(1),r2d(thisPolarAngle(1)),thisEccentricity(1),thisRfHalfWidth(1),a.params.pRFFit.rfType,thisx(1),thisy(1),num2str(r(1),'%0.2f '), x,y,z,thisR2(2),r2d(thisPolarAngle(2)),thisEccentricity(2),thisRfHalfWidth(2),a.params.pRFFit.rfType,thisx(2),thisy(2),num2str(r(2),'%0.2f ')));
% plot the rf
for scanNum = 1:2
a = subplot(6,5,15+5*(scanNum-1));
imagesc(d{scanNum}.stimX(:,1),d{scanNum}.stimY(1,:),flipud(m{scanNum}.rfModel'));
set(a,'Box','off');
set(a,'Color',[0.8 0.8 0.8]);
set(a,'TickDir','out');
axis equal
axis tight
hold on
hline(0,'w:');vline(0,'w:');
if scanNum ==1
    title('Scan 1: BarsTask')
else
    title ('Scan 2: BarsTaskFixation')
end
end
% plot the canonical
subplot(5,5,5);cla
plot(m{1}.canonical.time,m{1}.canonical.hrf,'k-');
title(sprintf('lag: %0.2f tau: %0.2f',m{1}.p.canonical.timelag,m{1}.p.canonical.tau));

% display the stimulus images
plotStim(gpRFPlot.t);

% now set callback
set(gpRFPlot.fignum,'WindowButtonMotionFcn',@pRFPlotMoveMouse);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrprFPloMoveMouse    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pRFPlotMoveMouse(hWindow,event)

global gpRFPlot;
if ~ishandle(gpRFPlot.a),return,end

currentPoint = get(gpRFPlot.a ,'CurrentPoint');
coord = round(currentPoint(1,1:2));
a = axis(gpRFPlot.a);
if (coord(1) >= a(1)) && (coord(1) <= a(2)) && (coord(2) >= a(3)) && (coord(2) <= a(4))
  % move rectangle
  pos = get(gpRFPlot.hRect,'Position');
  pos(1) = coord(1)-4;
  set(gpRFPlot.hRect,'Position',pos);
  % redisplay stimulus images
  plotStim(coord(1))
end

%%%%%%%%%%%%%%%%%%
%    plotStim    %
%%%%%%%%%%%%%%%%%%
function plotStim(t)

global gpRFPlot;

for scanNum = 1:2
for i = 1:5
  a = subplot(6,5,20+5*(scanNum-1)+i,'Parent',gpRFPlot.fignum);
  cla(a);
  thist = t-5+i;
  if thist >= 1 
    im = [];
    % get the scan and volume
    thisScan = gpRFPlot.d{scanNum}.concatInfo.whichScan(thist);
    thisVolume = gpRFPlot.d{scanNum}.concatInfo.whichVolume(thist);
    junkFrames = gpRFPlot.d{scanNum}.concatInfo.totalJunkedFrames(thisScan);
    im(:,:,3) = flipud(0.7*gpRFPlot.d{scanNum}.stim{thisScan}.im(:,:,thisVolume+junkFrames)');
    im(:,:,2) = flipud(0.7*gpRFPlot.d{scanNum}.stim{thisScan}.im(:,:,thisVolume+junkFrames)');
    im(:,:,1) = flipud(0.7*gpRFPlot.d{scanNum}.stim{thisScan}.im(:,:,thisVolume+junkFrames)'+0.3*gpRFPlot.rfModel{scanNum}');
    % swap and flip so that it will display correctly
    image(gpRFPlot.d{scanNum}.stimX(:,1),gpRFPlot.d{scanNum}.stimY(1,:),im,'Parent',a);
    axis image
    hold(a,'on');
    hline(0,'w:',a);vline(0,'w:',a);
    title(a,sprintf('t=%i',thist));
  end
end
end

%%%%%%%%%%%%%%%%%%%%
%    pRFPlotROI    %
%%%%%%%%%%%%%%%%%%%%
function pRFPlotROI(v,roi,d,a,r2,eccentricity,polarAngle,rfHalfWidth)
for scanNum = 1:2
if length(roi)
  % check for already plotted
  minr2 = viewGet(v,'overlayMin','r2');
  %scanNum = viewGet(v,'curScan');
  groupNum = viewGet(v,'curGroup');
  global gpRFPlotROI
  checkParams = {'roi','minr2','a','scanNum','groupNum'};
  replot = false;
  % if shift key is down then replot
  f = viewGet(v,'fignum');
 
  if ~isempty(f) && any(strcmp(get(f,'CurrentModifier'),'shift')),replot=true;end
  for i = 1:length(checkParams)
    if ~isfield(gpRFPlotROI,checkParams{i}) || ~isequal(gpRFPlotROI.(checkParams{i}),eval(checkParams{i}))
      replot = true;
    end
    gpRFPlotROI.(checkParams{i}) = eval(checkParams{i});
  end
  if ~replot, return, end
  disp(sprintf('(pRFPlot) Displaying ROI fig'));
  mlrSmartfig('pRFPlotROI','reuse');%clf
  
  minX(scanNum) = min(d{scanNum}.stimX(:));
  maxX(scanNum) = max(d{scanNum}.stimX(:));
  minY(scanNum) = min(d{scanNum}.stimY(:));
  maxY(scanNum) = max(d{scanNum}.stimY(:));
  
  % see what kind of fit we have.
  if strcmp(a.params.pRFFit.rfType,'gaussian-hdr')
    % plot also the hdr parameters
    numRowsPerROI = 2;
    numCols = 3;
    % set up fields for plotting extra hdr parameters
    if a.params.pRFFit.diffOfGamma
      plotParams = [4 5 6 7 8];
      plotParamsNames = {'timelag','tau','amplitudeRatio','timelag2','tau2'};
      numCols = 5;
    else
      plotParams = [4 5];
      plotParamsNames = {'timelag','tau'};
    end
  else
    numRowsPerROI = 1;
    numCols = 3;
    plotParams = [];
    plotParamsNames = {};
  end
 
  
  for roiNum = 1:length(roi)
    % get coordinates
    % roiCoords = getROICoordinates(v,roi{roiNum},[],[],'straightXform=1');
    roiCoords{scanNum}{roiNum} = getROICoordinates(v,roi{roiNum});
    roiCoordsLinear{scanNum}{roiNum} = sub2ind(viewGet(v,'scanDims'),roiCoords{scanNum}{roiNum}(1,:),roiCoords{scanNum}{roiNum}(2,:),roiCoords{scanNum}{roiNum}(3,:));
    % get values for the roi
    thisr2{scanNum}{roiNum} = r2{scanNum}(roiCoordsLinear{scanNum}{roiNum});
    % only use voxels above current r2 min
    roiCoordsLinear{scanNum}{roiNum} = roiCoordsLinear{scanNum}{roiNum}(find(thisr2{scanNum}{roiNum} >minr2));
    % sort them
    [thisr2sorted{scanNum}{roiNum} r2index{scanNum}{roiNum}] = sort(r2{scanNum}(roiCoordsLinear{scanNum}{roiNum}));
    roiCoordsLinear{scanNum}{roiNum} = roiCoordsLinear{scanNum}{roiNum}(r2index{scanNum}{roiNum});
    % get values for these voxels
    thisr2{scanNum}{roiNum} = r2{scanNum}(roiCoordsLinear{scanNum}{roiNum});
    thisEccentricity{scanNum}{roiNum} = eccentricity{scanNum}(roiCoordsLinear{scanNum}{roiNum});
    thisPolarAngle{scanNum}{roiNum} = polarAngle{scanNum}(roiCoordsLinear{scanNum}{roiNum});
    thisRfHalfWidth{scanNum}{roiNum} = rfHalfWidth{scanNum}(roiCoordsLinear{scanNum}{roiNum});
    % convert to cartesian
    [thisX{scanNum}{roiNum} thisY{scanNum}{roiNum}] = pol2cart(thisPolarAngle{scanNum}{roiNum},thisEccentricity{scanNum}{roiNum});
    c = [1 1 1];
    
    % plot RF coverage
    subplot(2*length(roi)*numRowsPerROI,numCols,1+(scanNum*roiNum-1)*numCols*numRowsPerROI);
    for i = 1:length(thisX{scanNum}{roiNum})
      if ~isnan(thisr2{scanNum}{roiNum}(i))
        plotCircle(thisX{scanNum}{roiNum}(i),thisY{scanNum}{roiNum}(i),thisRfHalfWidth{scanNum}{roiNum}(i),1-c*thisr2{scanNum}{roiNum}(i)/max(thisr2{scanNum}{roiNum}));
        hold on
      end
    end
    xaxis(minX(scanNum),maxX(scanNum));
    yaxis(minY(scanNum),maxY(scanNum));
    axis square
    hline(0);
    vline(0);
    xlabel('x (deg)');
    ylabel('y (deg)');
    title(sprintf('Scan %d: %s rf (r2 cutoff: %0.2f)',scanNum,roi{roiNum}.name,minr2));
    
    % plot RF centers
    subplot(2*length(roi)*numRowsPerROI,numCols,2+(scanNum*roiNum-1)*numCols*numRowsPerROI);
    for i = 1:length(thisX{scanNum}{roiNum})
      if ~isnan(thisr2{scanNum}{roiNum}(i))
        plot(thisX{scanNum}{roiNum}(i),thisY{scanNum}{roiNum}(i),'k.','Color',1-c*thisr2{scanNum}{roiNum}(i)/max(thisr2{scanNum}{roiNum}), 'markersize', 10);
        hold on
      end
    end
    xaxis(minX(scanNum),maxX(scanNum));
    yaxis(minY(scanNum),maxY(scanNum));
    axis square
    hline(0);
    vline(0);
    xlabel('x (deg)');
    ylabel('y (deg)');
    title(sprintf('Scan %d: %s centers',scanNum,roi{roiNum}.name));
    
    % plot eccentricity vs. rfHalfWidth
    subplot(2*length(roi)*numRowsPerROI,numCols,3+(scanNum*roiNum-1)*numCols*numRowsPerROI);
    for i = 1:length(thisX{scanNum}{roiNum})
      if ~isnan(thisr2{scanNum}{roiNum}(i))
        plot(thisEccentricity{scanNum}{roiNum}(i),thisRfHalfWidth{scanNum}{roiNum}(i),'k.','Color',1-c*thisr2{scanNum}{roiNum}(i)/max(thisr2{scanNum}{roiNum}), 'markersize', 10);
        hold on
      end
    end
    hold on
    % limit the fit to the central 6 deg (b/c it is often off for higher eccentricities)
    eccLimit = 6;
    ind = thisEccentricity{scanNum}{roiNum} <= eccLimit;
    if any(ind)
%      regfit = myregress(thisEccentricity(ind),thisRfHalfWidth(ind),0,0);
      w = diag(thisr2{scanNum}{roiNum}(ind));
      x = thisEccentricity{scanNum}{roiNum}(ind);
      x = [x(:) ones(size(x(:)))];
      y = thisRfHalfWidth{scanNum}{roiNum}(ind);
      beta = ((x'*w*x)^-1)*(x'*w)*y';
      maxXaxis = min(maxX(scanNum),maxY(scanNum));
      xaxis(0,maxXaxis);
      yaxis(0,maxXaxis);
      if ~isempty(beta)
	plot([0 maxXaxis],[0 maxXaxis]*beta(1)+beta(2),'k-');
      end
      xlabel('Eccentricity (deg)');
      ylabel('RF half width (deg)');
%      title(sprintf('slope: %0.2f (%s) offset: %0.2f (%s) (r2=%0.2f)',beta(1),pvaldisp(regfit.pm),beta(2),pvaldisp(regfit.pb),regfit.r2));
      axis square
    else
      disp(sprintf('(pRFPlot) No matching fits to plot with eccentricity less than %f',eccLimit));
    end
    title(sprintf('Scan %d: %s', scanNum,roi{roiNum}.name))
    % plot hdr parameters, first get the voxels to plot
    [temp{scanNum} dCoords{scanNum}] = intersect(d{scanNum}.linearCoords,roiCoordsLinear{scanNum}{roiNum});
    for i = 1:length(plotParams)
      subplot(2*length(roi)*numRowsPerROI,numCols,numCols+i+(scanNum*roiNum-1)*numCols*numRowsPerROI);
      hist(d{scanNum}.params(plotParams(i),dCoords{scanNum}));
      xlabel(plotParamsNames{i});
      ylabel('n');
      if exist('plotmean')==2
        plotmean(d{scanNum}.params(plotParams(i),dCoords{scanNum}));
      end
    end
  end
end
end


%%%%%%%%%%%%%%%%%%%%
%    plotCircle    %
%%%%%%%%%%%%%%%%%%%%
function h = plotCircle(xCenter,yCenter,radius,c)

a = 0:0.01:2*pi;
h = plot(xCenter+radius*cos(a),yCenter+radius*sin(a),'k-','Color',c);


%%%%%%%%%%%%%
%%   r2d   %%
%%%%%%%%%%%%%
function degrees = r2d(angle)

degrees = (angle/(2*pi))*360;

% if larger than 360 degrees then subtract
% 360 degrees
while (sum(degrees>360))
  degrees = degrees - (degrees>360)*360;
end

% if less than 360 degreees then add 
% 360 degrees
while (sum(degrees<-360))
  degrees = degrees + (degrees<-360)*360;
end

%% pRFCompareROI
%%%%%%%%%%%%%%%%%%%%
function pRFCompareROI(v,roi,d,a,r2,eccentricity,polarAngle,rfHalfWidth)
if length(roi)
global pRFPlotCompareROI
disp(sprintf('(pRFPlot) Displaying ROI fig: comparison'));
mlrSmartfig('pRFPlotCompareROI','reuse'); clf;
nScans = length(d);
ncols = 5;
for scanNum = 1:nScans
    minX(scanNum) = min(d{scanNum}.stimX(:));
    maxX(scanNum) = max(d{scanNum}.stimX(:));
    minY(scanNum) = min(d{scanNum}.stimY(:));
    maxY(scanNum) = max(d{scanNum}.stimY(:));

    roiCoords{scanNum}= getROICoordinates(v,roi{1});
    roiCoordsLinear{scanNum} = sub2ind(viewGet(v,'scanDims'),roiCoords{scanNum}(1,:),roiCoords{scanNum}(2,:),roiCoords{scanNum}(3,:));
    % get values for the roi
    thisr2{scanNum} = r2{scanNum}(roiCoordsLinear{scanNum});
    
    % get values for these voxels
    thisr2{scanNum} = r2{scanNum}(roiCoordsLinear{scanNum});
    thisEccentricity{scanNum} = eccentricity{scanNum}(roiCoordsLinear{scanNum});
    thisPolarAngle{scanNum} = polarAngle{scanNum}(roiCoordsLinear{scanNum});
    thisRfHalfWidth{scanNum} = rfHalfWidth{scanNum}(roiCoordsLinear{scanNum});
    % convert to cartesian
    [thisX{scanNum} thisY{scanNum}] = pol2cart(thisPolarAngle{scanNum},thisEccentricity{scanNum});
    c = [1 1 1];
end
numROI = length(roi);
 % plot r^2
    subplot(numROI,ncols,1+(numROI-1)*ncols);
    for i = 1:length(thisr2{2})
      if ~isnan(thisr2{2}(i)) && ~isnan(thisr2{1}(i))
        plot(thisr2{2}(i),thisr2{1}(i),'k.','Color',1-c*mean([thisr2{1}(i) thisr2{2}(i)])/max(mean([thisr2{1}; thisr2{2}])), 'markersize', 10);
        hold on
      end
    end
%     plot(thisr2{2}, thisr2{1}, 'k.', 'markersize', 10); hold on
%     lsline;
    x = 0:.005:1; y=x;
    plot(x,y,':', 'Color', [0.4,0.4,0.4])
    maxr2 = max([thisr2{2} thisr2{1}]);
    xaxis(0, round(maxr2*10)/10);
    yaxis(0, round(maxr2*10)/10);
    axis square
    xlabel('r^2: (Ave2)BarsTaskFixation')
    ylabel('r^2: (Ave1)BarsTask')
    title(sprintf('%s r^2',roi{1}.name));
        
 % plot RF center (x)
    subplot(numROI,ncols,2+(numROI-1)*ncols)
    for i = 1:length(thisX{2})
      if ~isnan(thisr2{2}(i)) && ~isnan(thisr2{1}(i))
        plot(thisX{2}(i),thisX{1}(i),'k.','Color',1-c*mean([thisr2{1}(i) thisr2{2}(i)])/max(mean([thisr2{1}; thisr2{2}])), 'markersize', 10);
        hold on
      end
    end
%     plot(thisX{2}, thisX{1}, 'k.', 'markersize', 10); hold on
%     lsline;
    x = minX(2):.5:maxX(2); y=x;
    plot(x,y,':', 'Color', [0.4,0.4,0.4])
    xaxis(min(minX(2),minX(1)), max(maxX(1),maxX(2)));
    yaxis(min(minX(1),minX(2)), max(maxX(1),maxX(2)));
    axis square
    xlabel('x: (Ave2)BarsTaskFixation')
    ylabel('x: (Ave1)BarsTask')
    title(sprintf('%s RF center (x)',roi{1}.name));
    
 % plot RF center (y)
    subplot(numROI,ncols,3+(numROI-1)*ncols)
    for i = 1:length(thisY{2})
      if ~isnan(thisr2{2}(i)) && ~isnan(thisr2{1}(i))
        plot(thisY{2}(i),thisY{1}(i),'k.','Color',1-c*mean([thisr2{1}(i) thisr2{2}(i)])/max(mean([thisr2{1}; thisr2{2}])), 'markersize', 10);
        hold on
      end
    end
%     plot(thisY{2}, thisY{1}, 'k.', 'markersize', 10); hold on
%     lsline;
    x = minY(2):.5:maxY(2); y=x;
    plot(x,y,':', 'Color', [0.4,0.4,0.4])
    xaxis(min(minY(2),minY(1)), max(maxY(1),maxY(2)));
    yaxis(min(minY(1),minY(2)), max(maxY(1),maxY(2)));
    axis square
    xlabel('y: (Ave2)BarsTaskFixation')
    ylabel('y: (Ave1)BarsTask')
    title(sprintf('%s RF center (y)',roi{1}.name));
    
 % plot rfHalfWidth
    subplot(numROI,ncols,4+(numROI-1)*ncols)
    for i = 1:length(thisRfHalfWidth{2})
      if ~isnan(thisr2{2}(i)) && ~isnan(thisr2{1}(i))
        plot(thisRfHalfWidth{2}(i),thisRfHalfWidth{1}(i),'k.','Color',1-c*mean([thisr2{1}(i) thisr2{2}(i)])/max(mean([thisr2{1}; thisr2{2}])), 'markersize', 10);
        hold on
      end
    end
%     plot(thisRfHalfWidth{2}, thisRfHalfWidth{1}, 'k.', 'markersize', 10); hold on
%     lsline
    x = min([thisRfHalfWidth{2} thisRfHalfWidth{1}]):.5:max([thisRfHalfWidth{2} thisRfHalfWidth{1}]); y=x;
    plot(x,y,':', 'Color', [0.4,0.4,0.4])
     %xaxis(min([thisRfHalfWidth{1} thisRfHalfWidth{2}]), max([thisRfHalfWidth{1} thisRfHalfWidth{2}]))
     %yaxis(min([thisRfHalfWidth{1} thisRfHalfWidth{2}]), max([thisRfHalfWidth{1} thisRfHalfWidth{2}]))
     xaxis(0,20)
     yaxis(0,20)
    axis square
    xlabel('rfHalfWidth: (Ave2)BarsTaskFixation')
    ylabel('rfHalfWidth: (Ave1)BarsTask')
    title(sprintf('%s RF half width',roi{1}.name));
    
 % plot eccentricity
    subplot(numROI,ncols,5+(numROI-1)*ncols)
    for i = 1:length(thisEccentricity{2})
      if ~isnan(thisr2{2}(i)) && ~isnan(thisr2{1}(i))
        plot(thisEccentricity{2}(i),thisEccentricity{1}(i),'k.','Color',1-c*mean([thisr2{1}(i) thisr2{2}(i)])/max(mean([thisr2{1}; thisr2{2}])), 'markersize', 10);
        hold on
      end
    end
%     plot(thisEccentricity{2}, thisEccentricity{1}, 'k.', 'markersize', 10); hold on
%     lsline
    x = min([thisEccentricity{2} thisEccentricity{1}]):1:max([thisEccentricity{2} thisEccentricity{1}]); y=x;
    plot(x,y,':', 'Color', [0.4,0.4,0.4])
%     xaxis(min([thisEccentricity{1} thisEccentricity{2}]), max([thisEccentricity{1} thisEccentricity{2}]))
%     yaxis(min([thisEccentricity{1} thisEccentricity{2}]), max([thisEccentricity{1} thisEccentricity{2}]))
     xaxis(0,50)
     yaxis(0,50)
    axis square
    xlabel('Eccentricity: (Ave2)BarsTaskFixation')
    ylabel('Eccentricity: (Ave1)BarsTask')
    title(sprintf('%s Eccentricity',roi{1}.name));
    

end

%% pRFCompareROI2
%%%%%%%%%%%%%%%%%%%%
function pRFCompareROI2(v,roi,d,a,r2,eccentricity,polarAngle,rfHalfWidth)
if length(roi)
global pRFPlotCompareROI2
disp(sprintf('(pRFPlot) Displaying ROI fig: comparison'));
mlrSmartfig('pRFPlotCompareROI2','reuse'); clf;
nScans = length(d);
ncols = 5;
for scanNum = 1:nScans
    minX(scanNum) = min(d{scanNum}.stimX(:));
    maxX(scanNum) = max(d{scanNum}.stimX(:));
    minY(scanNum) = min(d{scanNum}.stimY(:));
    maxY(scanNum) = max(d{scanNum}.stimY(:));

    roiCoords{scanNum}= getROICoordinates(v,roi{1});
    roiCoordsLinear{scanNum} = sub2ind(viewGet(v,'scanDims'),roiCoords{scanNum}(1,:),roiCoords{scanNum}(2,:),roiCoords{scanNum}(3,:));
    % get values for the roi
    thisr2{scanNum} = r2{scanNum}(roiCoordsLinear{scanNum});
    
    % get values for these voxels
    thisr2{scanNum} = r2{scanNum}(roiCoordsLinear{scanNum});
    thisEccentricity{scanNum} = eccentricity{scanNum}(roiCoordsLinear{scanNum});
    thisPolarAngle{scanNum} = polarAngle{scanNum}(roiCoordsLinear{scanNum});
    thisRfHalfWidth{scanNum} = rfHalfWidth{scanNum}(roiCoordsLinear{scanNum});
    % convert to cartesian
    [thisX{scanNum} thisY{scanNum}] = pol2cart(thisPolarAngle{scanNum},thisEccentricity{scanNum});
    c = [1 1 1];
end
numROI = length(roi);
 % plot r^2
  for rr = 1:2
    subplot(2*numROI,ncols,1+(numROI-1+rr-1)*ncols);
    for i = 1:length(thisr2{2})
      if ~isnan(thisr2{2}(i)) && ~isnan(thisr2{1}(i))
        plot(thisr2{2}(i),thisr2{1}(i),'k.','Color',1-c*thisr2{rr}(i)/max(thisr2{rr}), 'markersize', 10);
        hold on
      end
    end
%     plot(thisr2{2}, thisr2{1}, 'k.', 'markersize', 10); hold on
%     lsline;
    x = 0:.005:1; y=x;
    plot(x,y,':', 'Color', [0.4,0.4,0.4])
    maxr2 = max([thisr2{2} thisr2{1}]);
    xaxis(0, round(maxr2*10)/10);
    yaxis(0, round(maxr2*10)/10);
    axis square
    xlabel('r^2: (Ave2)BarsTaskFixation')
    ylabel('r^2: (Ave1)BarsTask')
    title(sprintf('%s r^2',roi{1}.name));
  end
        
 % plot RF center (x)
  for rr = 1:2
    subplot(2*numROI,ncols,2+(numROI-1+rr-1)*ncols)
    for i = 1:length(thisX{2})
      if ~isnan(thisr2{2}(i)) && ~isnan(thisr2{1}(i))
        plot(thisX{2}(i),thisX{1}(i),'k.','Color',1-c*thisr2{rr}(i)/max(thisr2{rr}), 'markersize', 10);
        hold on
      end
    end
%     plot(thisX{2}, thisX{1}, 'k.', 'markersize', 10); hold on
%     lsline;
    x = minX(2):.5:maxX(2); y=x;
    plot(x,y,':', 'Color', [0.4,0.4,0.4])
    xaxis(min(minX(2),minX(1)), max(maxX(1),maxX(2)));
    yaxis(min(minX(1),minX(2)), max(maxX(1),maxX(2)));
    axis square
    xlabel('x: (Ave2)BarsTaskFixation')
    ylabel('x: (Ave1)BarsTask')
    title(sprintf('%s RF center (x)',roi{1}.name));
  end
    
 % plot RF center (y)
  for rr = 1:2
    subplot(2*numROI,ncols,3+(numROI-1+rr-1)*ncols)
    for i = 1:length(thisY{2})
      if ~isnan(thisr2{2}(i)) && ~isnan(thisr2{1}(i))
        plot(thisY{2}(i),thisY{1}(i),'k.','Color',1-c*thisr2{rr}(i)/max(thisr2{rr}), 'markersize', 10);
        hold on
      end
    end
%     plot(thisY{2}, thisY{1}, 'k.', 'markersize', 10); hold on
%     lsline;
    x = minY(2):.5:maxY(2); y=x;
    plot(x,y,':', 'Color', [0.4,0.4,0.4])
    xaxis(min(minY(2),minY(1)), max(maxY(1),maxY(2)));
    yaxis(min(minY(1),minY(2)), max(maxY(1),maxY(2)));
    axis square
    xlabel('y: (Ave2)BarsTaskFixation')
    ylabel('y: (Ave1)BarsTask')
    title(sprintf('%s RF center (y)',roi{1}.name));
  end
      
 % plot rfHalfWidth
 for rr = 1:2
    subplot(2*numROI,ncols,4+(numROI-1+rr-1)*ncols)
    for i = 1:length(thisRfHalfWidth{2})
      if ~isnan(thisr2{2}(i)) && ~isnan(thisr2{1}(i))
        plot(thisRfHalfWidth{2}(i),thisRfHalfWidth{1}(i),'k.','Color',1-c*thisr2{rr}(i)/max(thisr2{rr}), 'markersize', 10);
        hold on
      end
    end
%     plot(thisRfHalfWidth{2}, thisRfHalfWidth{1}, 'k.', 'markersize', 10); hold on
%     lsline
    x = min([thisRfHalfWidth{2} thisRfHalfWidth{1}]):.5:max([thisRfHalfWidth{2} thisRfHalfWidth{1}]); y=x;
    plot(x,y,':', 'Color', [0.4,0.4,0.4])
     %xaxis(min([thisRfHalfWidth{1} thisRfHalfWidth{2}]), max([thisRfHalfWidth{1} thisRfHalfWidth{2}]))
     %yaxis(min([thisRfHalfWidth{1} thisRfHalfWidth{2}]), max([thisRfHalfWidth{1} thisRfHalfWidth{2}]))
     xaxis(0,20)
     yaxis(0,20)
    axis square
    xlabel('rfHalfWidth: (Ave2)BarsTaskFixation')
    ylabel('rfHalfWidth: (Ave1)BarsTask')
    title(sprintf('%s RF half width',roi{1}.name));
 end
    
 % plot eccentricity
  for rr = 1:2
    subplot(2*numROI,ncols,5+(numROI-1+rr-1)*ncols)
    for i = 1:length(thisEccentricity{2})
      if ~isnan(thisr2{2}(i)) && ~isnan(thisr2{1}(i))
        plot(thisEccentricity{2}(i),thisEccentricity{1}(i),'k.','Color',1-c*thisr2{rr}(i)/max(thisr2{rr}), 'markersize', 10);
        hold on
      end
    end
%     plot(thisEccentricity{2}, thisEccentricity{1}, 'k.', 'markersize', 10); hold on
%     lsline
    x = min([thisEccentricity{2} thisEccentricity{1}]):1:max([thisEccentricity{2} thisEccentricity{1}]); y=x;
    plot(x,y,':', 'Color', [0.4,0.4,0.4])
%     xaxis(min([thisEccentricity{1} thisEccentricity{2}]), max([thisEccentricity{1} thisEccentricity{2}]))
%     yaxis(min([thisEccentricity{1} thisEccentricity{2}]), max([thisEccentricity{1} thisEccentricity{2}]))
     xaxis(0,50)
     yaxis(0,50)
    axis square
    xlabel('Eccentricity: (Ave2)BarsTaskFixation')
    ylabel('Eccentricity: (Ave1)BarsTask')
    title(sprintf('%s Eccentricity',roi{1}.name));
  end

end

