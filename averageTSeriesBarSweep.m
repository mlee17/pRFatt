% averageTSeriesBarSweep.m
%
%      usage: 
%         by: Minyoung Lee
%       date: 3/26/16
%    purpose:
%
function averageTSeriesBarSweep(v, group, roiname, minEccentricity, maxEccentricity, maxpRFSize, minr2)

% which ROIs
% if ieNotDefined('roiname')
%     % ROIsFull = {'V1','V2v','V2d','V3v','V3d','V3A','V3B','LO1','LO2','V4','MT','IPS0','IPS1','IPS2','IPS3','IPS4','IPS5'};
%     roiname = {'V1','V2v','V2d','V3v','V3d','V3A','V3B','LO1','LO2','V4','IPS0','IPS1','IPS2','IPS3','IPS4'};
% end
hemi = {'l','r'};

roiname={'V1','V2v','V3A','IPS0'};
%%
if ~exist('mrParamsDialog')
  disp(sprintf('(---) You must have mrTools in your path to run this'));
  return
end

if ieNotDefined('minEccentricity') || ieNotDefined('maxEccentricity') || ieNotDefined('maxpRFSize') || ieNotDefined('minr2')
    disp(sprintf('Not enough input arguments'))
    return
end

% Specify view
if ieNotDefined('v')
    v = newView;
end

% Group and Scans
if ieNotDefined('group')
    group = 'Averages';
end    
v = viewSet(v, 'currentGroup', group);
nScans = viewGet(v, 'nScans', 'Averages');
if nScans ~= 2
    disp(sprintf('(---) There must be 2 scans in Averages'));
    return
end
barsTaskScanNum = viewGet(v,'scanNumFromDescription', 'BarsTask ');
if barsTaskScanNum ~= 1
    scanOrder = askuser(sprintf('ScanNum for barsTask: %d\nScan Info might be inaccurate. Do you wish to continue?', barsTaskScanNum),1);
    if ~scanOrder
        return;
    end
end

analysisFile = dir('Averages/pRFAnal/pRF_ave*');
v = loadAnalysis(v, ['pRFAnal/' analysisFile.name]);
%%
% get the d
d=cell(1,2);
for scanNum = 1:2
    d{scanNum} = viewGet(v,'d',scanNum);
    if isempty(d{scanNum}),disp(sprintf('(pRFPlot) Could not find d structure for this scan %d', scanNum));return,end
    % get the parametrs of the pRF fit
    r2{scanNum} = viewGet(v,'overlayData',scanNum,viewGet(v,'overlayNum','r2'));
    if isempty(r2{scanNum})
        disp(sprintf('(pRFPlot) pRF analysis has not been run on this scan %d', scanNum));
    return
    end
    polarAngle{scanNum} = viewGet(v,'overlayData',scanNum,viewGet(v,'overlayNum','polarAngle'));
    eccentricity{scanNum} = viewGet(v,'overlayData',scanNum,viewGet(v,'overlayNum','eccentricity'));
    rfHalfWidth{scanNum} = viewGet(v,'overlayData',scanNum,viewGet(v,'overlayNum','rfHalfWidth'));
    
    % get the params that have been run
    scanDims{scanNum} = viewGet(v,'scanDims',scanNum);
end
    subject = viewGet(v, 'subject');

% ROI file names
% prepare an empty cell array
roiFileNames = cell(length(hemi)*length(roiname),1);
for roi = 1:length(roiname)
    for h = 1:length(hemi)
        roiFileNames{2*roi+(h-2)} = [hemi{h} roiname{roi} '_C'];
    end
end
  
% Load ROIs and adds them to view
v = loadROI(v, roiFileNames);
% Load ROI tSeries
% tseries: cell array (nROIs x nScans), each element of which is nFrames x nVoxels matrix

tsROI = tseriesROI(v);

[barList, barTPoints] = getBarTPoints(v, d);


global fig180; global fig90; global fig270;
for roi = 1:length(roiname)

for h = 1:length(hemi)
    clearvars this* centeredTseries r2pass eccpass sizepass incld roiCoordsLinear roiCoords
    
    name = [hemi{h} roiname{roi} '_C'];

    roiCoords= getROICoordinates(v,name);
    roiCoordsLinear = sub2ind(scanDims{1},roiCoords(1,:),roiCoords(2,:),roiCoords(3,:));
    
    for scanNum = 1:nScans
           % get values for the roi
            thisr2{scanNum} = r2{scanNum}(roiCoordsLinear);
            thisEccentricity{scanNum} = eccentricity{scanNum}(roiCoordsLinear);
            thisRfHalfWidth{scanNum} = rfHalfWidth{scanNum}(roiCoordsLinear);
    end
    minEccentricity, maxEccentricity, maxpRFSize, minr2
    
    r2pass = find(thisr2{1}>=minr2 & thisr2{2}>=minr2);
    eccpass = intersect(find(thisEccentricity{1}>=minEccentricity & thisEccentricity{1}<=maxEccentricity), find(thisEccentricity{2}>=minEccentricity & thisEccentricity{2}<=maxEccentricity));
    %find(mean([thisEccentricity{1}; thisEccentricity{2}])<=maxEccentricity);
%     eccpassLower = intersect(find(thisEccentricity{1}>=minEccentricity), find(thisEccentricity{2}>=minEccentricity));
    %find(mean([thisEccentricity{1}; thisEccentricity{2}])>=minEccentricity);
    sizepass = find(thisRfHalfWidth{1}<maxpRFSize & thisRfHalfWidth{2}<maxpRFSize);
    incld = intersect(r2pass,intersect(eccpass, sizepass));
    
    roiCoordsLinear = roiCoordsLinear(incld);
    
    for scanNum = 1:nScans    
          
        % get values for these voxels
            thisr2{scanNum} = r2{scanNum}(roiCoordsLinear);
            thisEccentricity{scanNum} = eccentricity{scanNum}(roiCoordsLinear);
            thisPolarAngle{scanNum} = polarAngle{scanNum}(roiCoordsLinear);
            thisRfHalfWidth{scanNum} = rfHalfWidth{scanNum}(roiCoordsLinear);
        % convert to cartesian
            [thisX{scanNum} thisY{scanNum}] = pol2cart(thisPolarAngle{scanNum},thisEccentricity{scanNum});
            c = [1 1 1];
        
        % ALSO, tSeries for the ROI!!!
        % tseries: cell array (nROIs x nScans), each element of which is nFrames x nVoxels matrix
            tsROI{2*roi+(h-2),scanNum} = tsROI{2*roi+(h-2),scanNum}(:,incld);
            
    end

    % rf centers in baseScan(Fixation Cond)
    % Find the timepoint when the bar crossing the rf Center in
    % *** FIXATION *** Scan !!!
    % for EACH BAR SWEEP !!!!!!
%     nVoxels = length(roiCoordsLinear);
    baseScan = 2; %fixation
    
    aveTPlotTaskFix(d,baseScan,h, hemi,roi,roiname, tsROI,thisX,thisY,thisr2, barList, barTPoints, 180)
    aveTPlotTaskFix(d,baseScan,h, hemi,roi, roiname,tsROI,thisX,thisY,thisr2, barList, barTPoints,90)
    aveTPlotTaskFix(d,baseScan,h, hemi,roi, roiname,tsROI,thisX,thisY,thisr2, barList, barTPoints,270)

end
end
        deleteView(v)
        

       
%%
function aveTPlotTaskFix(d,baseScan,h, hemi, roi,roiname, tsROI,thisX,thisY,thisr2, barList, barTPoints,angle)

global fig180; global fig90; global fig270;
        % BAR ANGLE 1: 180 right -> left (Horizontal movement / X position only)
       stimStart = barTPoints(1,barList == angle);
       stimEnd = barTPoints(2,barList == angle);
%        im = d{2}.stim{1}.im(:,1,stimStart:stimEnd);
%        squeeze(d{2}.stim{1}.im(:,1,stimStart:stimEnd))
       
       xt = squeeze(d{2}.stim{1}.im(:,1,:));
       yt = squeeze(d{2}.stim{1}.im(1,:,:));
       
    if any(angle == [180 0])
       stim = d{2}.stimX(:,1);
       thisXorY = thisX;
    elseif any(angle==[90 270])
       stim = d{2}.stimY(1,:);
       thisXorY = thisY;
    else
        return
    end
      
       
       stepStamp = zeros(1,2);
       stepStamp(1) = stimStart;
       ind = 2;
       for s = stimStart+1:stimEnd
           if any(angle == [180 0])
            if find(xt(:,s),1,'first') ~= find(xt(:,s-1),1,'first')
               stepStamp(ind) = s;
               ind=ind+1;
            end
           elseif any(angle==[90 270])
            if find(yt(:,s),1,'first') ~= find(yt(:,s-1),1,'first')
               stepStamp(ind) = s;
               ind=ind+1;
            end
           end
       end
       
       ind = 1;
       for s = 1:length(stepStamp)
           if any(angle == [180 0])
               meanPos(s) = mean(stim(find(xt(:,stepStamp(s)))));
           elseif any(angle==[90 270])
               meanPos(s) = mean(stim(find(yt(:,stepStamp(s)))));
           end

           xrangeUp = meanPos(s)+1.5;%stimX(find(xt(:,stepStamp(s)),1,'last')); % larger value
           xrangeLow = meanPos(s)-1.5;%stimX(find(xt(:,stepStamp(s)),1,'first')); % smaller value
           centerSweep{s} = find(thisXorY{baseScan} <= xrangeUp & thisXorY{baseScan} >= xrangeLow);
           for nvox = 1:length(centerSweep{s})
               centeredTseries.task.tseries{ind} = tsROI{2*roi+(h-2),1}(stepStamp(s)-30:stepStamp(s)+30,centerSweep{s}(nvox));
               centeredTseries.task.r2(ind) = thisr2{1}(centerSweep{s}(nvox));
               centeredTseries.task.weightedTseries{ind} = centeredTseries.task.tseries{ind} * centeredTseries.task.r2(ind);
               
               centeredTseries.fix.tseries{ind} = tsROI{2*roi+(h-2),2}(stepStamp(s)-30:stepStamp(s)+30,centerSweep{s}(nvox));
               centeredTseries.fix.r2(ind) = thisr2{2}(centerSweep{s}(nvox));
               centeredTseries.fix.weightedTseries{ind} = centeredTseries.fix.tseries{ind} * centeredTseries.fix.r2(ind);
               
               centeredTseries.task.weightedTseriesFixW{ind} = centeredTseries.task.tseries{ind} * centeredTseries.fix.r2(ind);
               ind = ind+1;
           end
       end
       
       % Weighted mean -> divide by the sum of the weights (r^2)
        meanTseries_task = sum([centeredTseries.task.weightedTseries{:}],2) / sum(centeredTseries.task.r2);
        meanTseries_fix = sum([centeredTseries.fix.weightedTseries{:}],2) / sum(centeredTseries.fix.r2);
        
        name = [hemi{h} roiname{roi}];
        
%         meanTseries_task_FW = sum([centeredTseries.task.weightedTseriesFixW{:}],2)/ sum(centeredTseries.fix.r2);
        if angle==180 ;mlrSmartfig('fig180', 'reuse'); elseif angle==90; mlrSmartfig('fig90', 'reuse'); elseif angle==270; mlrSmartfig('fig270', 'reuse'); end
        t = -15:.5:15;
        subplot(length(roiname),length(hemi),h+(length(hemi)*(roi-1)))
        plot(t, meanTseries_fix); hold on;
        plot(t, meanTseries_task, 'r-');
        xaxis(-15,15);
%         yaxis(-1,1);
        title(sprintf('%s %d', name, angle))
        


    

%%%%%%%%%%%%%%%%%%%%%%%%%
%%    getBarTPoints    %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function [barList, barTPoints] = getBarTPoints(v, d)

stimfile = viewGet(v,'stimfile', 2);
barAngle = stimfile{1}.task{2}{1}.block.parameter.barAngle;
%start and end times (1st & 2nd rows)
timepoints = zeros(2,length(barAngle));
volPerSec = 2; lenBlank = 10; lenBar = 24;
nVolsBlank = volPerSec*lenBlank; 
nVolsBar = volPerSec*lenBar;
% initial Blank
timepoints(1,1) = 1;
timepoints(2,1) = nVolsBlank;
for seq = 2:length(barAngle)
    timepoints(1,seq) = timepoints(2,seq-1) + 1;
    if barAngle(seq) == -1
        timepoints(2,seq) = timepoints(2,seq-1)+nVolsBlank;
    else
        timepoints(2,seq) = timepoints(2,seq-1)+nVolsBar;
    end
end
%this line may be fixed later for the correct number of volumes!!!!!
timepoints(2,length(barAngle)) = length(d{1}.stimT);

barList = barAngle(barAngle~=-1);
barTPoints = timepoints(:,barAngle~=-1);

%%%%%% barAngle %%%%% 

%   <cardinal>
%   0: left -> right
%   180: right -> left
%   90: downward
%   270: upward
%
%   <diagonal>
%   45: upper left -> lower right
%   135: upper right -> lower left
%   225: lower right -> upper left
%   315: lower left -> upper right

%%%%%%%%%%%%%%%%%%%%%








