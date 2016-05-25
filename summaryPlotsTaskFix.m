% summaryPlotsTaskFix.m
%
%    usage: ROIStruct = summaryPlotsTaskFix(v, varargin)
%       by: Minyoung Lee
%     date:
%  purpose:


% [ROIStruct, bothHemi] = summaryPlotsTaskFix([],0.2, 25,'C','collapse')


function [ROIStruct, bothHemi]= summaryPlotsTaskFix(v, r2cutoff,eccCutoff, ROItype, plottype, varargin)

if ~exist('mrParamsDialog')
  disp(sprintf('(---) You must have mrTools in your path to run this'));
  return
end

if ieNotDefined('v')
    v = newView;
end

% groupnum =viewGet(v, 'groupnum','Averages');
v = viewSet(v, 'currentGroup', 'Averages'); % "Average" Group

numScans = viewGet(v, 'nScans', 'Averages');
if numScans ~= 2
    disp(sprintf('(---) There must be 2 scans in Averages'));
    return
end

barsTaskScanNum = viewGet(v,'scanNumFromDescription', 'BarsTask ');
% barsTaskFixationScanNum = viewGet(v,'scanNumFromDescription', 'BarsTaskFixation');

if barsTaskScanNum ~= 1 %|| barsTaskFixationScanNum ~= 2 %isempty(barsTaskScanNum) || isempty(barsTaskFixationScanNum)
    scanOrder = askuser(sprintf('ScanNum for barsTask: %d\nScanNum for barsTaskFixation: %d\nScan Info might be inaccurate. Do you wish to continue?', barsTaskScanNum, barsTaskFixationScanNum),1);
    if ~scanOrder
        return;
    end
end

analysisFile = dir('Averages/pRFAnal/pRF_ave*');
v = loadAnalysis(v, ['pRFAnal/' analysisFile.name]);

hemi = {'l','r'};
%ROIs = {'V1','V2v','V2d','V3v','V3d','V3A','V3B','LO1','LO2','V4','MT','IPS0','IPS1','IPS2','IPS3','IPS4','IPS5'};
ROIs = {'V1','V2v','V2d','V3v','V3d','V3A','V3B','LO1','LO2','V4','IPS0','IPS1','IPS2','IPS3','IPS4'};%,'IPS5'};


% roiDir = viewGet(v,'roidir')
% stimfile = viewGet(v,'stimfile',1,3);

% % check if pRF has been run
% a = viewGet(v,'Analysis');
% if ~isfield(a,'type') || ~strcmp(a.type,'pRFAnal')
%   disp(sprintf('(pRFPlot) pRF analysis has not been run on this scan'));
%   return
% end

d=cell(1,2); r2=cell(1,2); polarAngle=cell(1,2); eccentricity=cell(1,2); rfHalfWidth=cell(1,2); scanDims=cell(1,2); whichVoxel=cell(1,2); r=cell(1,2); 
params=cell(1,2); paramsInfo=cell(1,2); %m=cell(1,2); 
for scanNum = 1:2
    % get the d
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

% x,y,z ERROR!!!!
% whichVoxel(scanNum) = find(d{scanNum}.linearCoords == sub2ind(scanDims{scanNum},x,y,z));
% r(scanNum) = d{scanNum}.r(whichVoxel(scanNum),:);
% 
% params{scanNum} = d{scanNum}.params(:,whichVoxel(scanNum));
% if isfield(d{scanNum},'paramsInfo')
%   paramsInfo{scanNum} = d{scanNum}.paramsInfo;
% else
%   paramsInfo{scanNum} = [];
% end
% 
% 
% % get params
% m{scanNum} = pRFFit(v,scanNum,x,y,z,'stim',d{scanNum}.stim,'getModelResponse=1','params',params{scanNum},'concatInfo',d{scanNum}.concatInfo,'fitTypeParams',a.params.pRFFit,'paramsInfo',paramsInfo{scanNum});

end

% roi
% if ~shiftDown
%   pRFPlotROI(v,roi,d,a,r2,eccentricity,polarAngle,rfHalfWidth);
% r2cutoff = 0.20;
[ROIStruct,bothHemi] = pRFPlotMultipleROIs(v,ROIs,hemi,r2cutoff, eccCutoff, ROItype,plottype,d,r2,eccentricity,polarAngle,rfHalfWidth);
%   pRFCompareROI2(v,roi,d,a,r2,eccentricity,polarAngle,rfHalfWidth);
% end

deleteView(v)



%%%%%%%%%%%%%%%%%%%%
%% pRFPlotMultipleROIs
%%%%%%%%%%%%%%%%%%%%
function [ROIStruct, bothHemi] = pRFPlotMultipleROIs(v,ROIs,hemi,r2cutoff, eccCutoff, ROItype,plottype,d,r2,eccentricity,polarAngle,rfHalfWidth)
 subject = viewGet(v, 'subject');
% roiGroupN = [5,6,6];
for h = 1:length(hemi)
nScans = length(d);
numROI = length(ROIs);
    
    for roi = 1:numROI
%         ROIStruct{roi}.
        roiName = sprintf('%s%s_%s',hemi{h}, ROIs{roi}, ROItype);
        clearvars roiCoords* this*
        roiCoords{1}= getROICoordinates(v,roiName);
        if isempty(roiCoords{1})
            continue
        end
        for scanNum = 1:nScans
%             minX(scanNum) = min(d{scanNum}.stimX(:));
%             maxX(scanNum) = max(d{scanNum}.stimX(:));
%             minY(scanNum) = min(d{scanNum}.stimY(:));
%             maxY(scanNum) = max(d{scanNum}.stimY(:));
%                     
            roiCoords{scanNum}= getROICoordinates(v,roiName);
            roiCoordsLinear{scanNum} = sub2ind(viewGet(v,'scanDims'),roiCoords{scanNum}(1,:),roiCoords{scanNum}(2,:),roiCoords{scanNum}(3,:));
        % get values for the roi
            orgthisr2{scanNum} = r2{scanNum}(roiCoordsLinear{scanNum});
            orgthisEccentricity{scanNum} = eccentricity{scanNum}(roiCoordsLinear{scanNum});
            orgthisRfHalfWidth{scanNum} = rfHalfWidth{scanNum}(roiCoordsLinear{scanNum});
        end
        
        if ieNotDefined('eccCutoff')
            eccCutoff = 25;
        end
        
        r2pass = intersect(find(orgthisr2{1}>=r2cutoff),find(orgthisr2{2}>=r2cutoff));
        eccpass = find(mean([orgthisEccentricity{1}; orgthisEccentricity{2}])<=eccCutoff);
        incld = intersect(r2pass, eccpass);
        if strcmp(ROIs{roi}, 'V1') || strcmp(ROIs{roi}, 'V2d') || strcmp(ROIs{roi}, 'V2v') || strcmp(ROIs{roi}, 'V3v') || strcmp(ROIs{roi}, 'V3d')
            %widthCutoff = 12.5;
            widthCutoff = 15;
            sizepass = find(orgthisRfHalfWidth{1} <= widthCutoff & orgthisRfHalfWidth{2} <= widthCutoff);%find(mean([orgthisRfHalfWidth{1} ; orgthisRfHalfWidth{2}]) <=widthCutoff);
            incld = intersect(incld, sizepass);
        end
            
        
        for scanNum=1:nScans
        % exclude voxels if their pRF models explained less than XX % of
        % response variance in either condition (only use voxels above the
        % cutoff value)
%             roiCoordsLinear{scanNum} = roiCoordsLinear{scanNum}(intersect(find(thisr2{1}>=r2cutoff),find(thisr2{2}>=r2cutoff)));
            roiCoordsLinear{scanNum} = roiCoordsLinear{scanNum}(incld);
        end
        
        for scanNum = 1:nScans    
          
        % get values for these voxels
            thisr2{scanNum} = r2{scanNum}(roiCoordsLinear{scanNum});
            thisEccentricity{scanNum} = eccentricity{scanNum}(roiCoordsLinear{scanNum});
            thisPolarAngle{scanNum} = polarAngle{scanNum}(roiCoordsLinear{scanNum});
            thisRfHalfWidth{scanNum} = rfHalfWidth{scanNum}(roiCoordsLinear{scanNum});
        % convert to cartesian
            [thisX{scanNum} thisY{scanNum}] = pol2cart(thisPolarAngle{scanNum},thisEccentricity{scanNum});
            
        end
           
    ROIStruct{h}{roi} = struct('Subject', subject, 'Scan', {'BarsTask','BarsTaskFixation'},'roiName', roiName, 'roiCoords', roiCoords, 'roiCoordsLinear', roiCoordsLinear, ...
        'thisr2', thisr2, 'thisX', thisX, 'thisY', thisY, 'thisEccentricity', thisEccentricity, ...
        'thisPolarAngle', thisPolarAngle, 'thisRfHalfWidth', thisRfHalfWidth, 'orgthisr2', orgthisr2, 'orgthisEccentricity', orgthisEccentricity, ...
        'orgthisRfHalfWidth', orgthisRfHalfWidth);
            
    end
end
    
ROIStruct{3}.r2cutoff = r2cutoff;
ROIStruct{3}.meanEccCutoff = eccCutoff;

        for hemi = 1:2
            for roi = 1:length(ROIs)
                if roi == length(ROIs) && length(ROIStruct{hemi}) < roi
                    continue
                end
            if isempty(ROIStruct{hemi}{roi})
            for type = 1:2
                ROIStruct{hemi}{roi}(type).roiName = [];
                ROIStruct{hemi}{roi}(type).thisr2 = []; ROIStruct{hemi}{roi}(type).thisX = []; ROIStruct{hemi}{roi}(type).thisY = []; 
                ROIStruct{hemi}{roi}(type).thisEccentricity = []; ROIStruct{hemi}{roi}(type).thisRfHalfWidth=[]; 
            end
            end
            end
        end
        

bothHemi = [];
    if strcmp(plottype,'separate') 
%         plotScatterROI(d,ROIs, hemi,'left', thisX, thisY, thisEccentricity, thisRfHalfWidth, r2cutoff)
%         plotScatterROI(d,ROIs, hemi,'right', thisr2, thisX, thisY, thisEccentricity, thisRfHalfWidth, r2cutoff)
    elseif strcmp(plottype, 'collapse')
        hemi = cell(1,2);
     
        for roi = 1:length(ROIs)
            for cond = 1:2
            bothHemi.thisr2{roi}{cond} = [ROIStruct{1}{roi}(cond).thisr2, ROIStruct{2}{roi}(cond).thisr2];
            bothHemi.thisX{roi}{cond} = [ROIStruct{1}{roi}(cond).thisX, ROIStruct{2}{roi}(cond).thisX];
            bothHemi.thisY{roi}{cond} = [ROIStruct{1}{roi}(cond).thisY, ROIStruct{2}{roi}(cond).thisY];
            bothHemi.thisEccentricity{roi}{cond} = [ROIStruct{1}{roi}(cond).thisEccentricity, ROIStruct{2}{roi}(cond).thisEccentricity];
            bothHemi.thisRfHalfWidth{roi}{cond} = [ROIStruct{1}{roi}(cond).thisRfHalfWidth, ROIStruct{2}{roi}(cond).thisRfHalfWidth];
            end
        end
        plotScatterROI(d,ROIs, hemi, 'both',bothHemi, r2cutoff)%bothHemi.thisr2, bothHemi.thisX, bothHemi.thisY, bothHemi.thisEccentricity, bothHemi.thisRfHalfWidth, r2cutoff)
    end


%%%%%%%%%%%%%%%%%%%%
%% plotScatterROI %%
%%%%%%%%%%%%%%%%%%%%
function plotScatterROI(d,ROIs,hemi,hemispec, datastruct, r2cutoff)%thisr2, thisX, thisY, thisEccentricity, thisRfHalfWidth, r2cutoff)
%global pRFPlotMultipleROIs;
disp(sprintf('(pRFPlot) Displaying ROI fig: comparing two conditions within each ROI'));
%mlrSmartfig('pRFPlotMultipleROIs','reuse'); clf;
ncols = 5;    
scrsz = get(0, 'ScreenSize'); 
c = [1 1 1];
for s = 1:2
    minX(s) = min(d{s}.stimX(:));
    maxX(s) = max(d{s}.stimX(:));
    minY(s) = min(d{s}.stimY(:));
    maxY(s) = max(d{s}.stimY(:));
end

if strcmp(hemispec,'right')
    h =2;
else
    h=1;
end

for roi = 1:length(ROIs)

    thisr2 = datastruct.thisr2{roi}; thisX = datastruct.thisX{roi}; thisY = datastruct.thisY{roi}; 
    thisEccentricity= datastruct.thisEccentricity{roi}; thisRfHalfWidth=datastruct.thisRfHalfWidth{roi};   

     
     
if roi == 1; figure('Position', [1 scrsz(4)/1.2 scrsz(3)/2 scrsz(4)/1.2]) ; clf; roiGroupN =5; subtract = roi; end
        if roi == 6; figure('Position', [1 scrsz(4)/1.2 scrsz(3)/2 scrsz(4)/1.2]); clf; roiGroupN =5; subtract = roi; end
        if roi == 11; figure('Position', [1 scrsz(4)/1.2 scrsz(3)/2 scrsz(4)/1.2]); clf; roiGroupN =6; subtract = roi; end    
        
        % plot r^2
        subplot(roiGroupN,ncols,1+(roi-subtract)*ncols);
        for i = 1:length(thisr2{2})
                plot(thisr2{2}(i),thisr2{1}(i),'k.','Color',1-c*mean([thisr2{1}(i) thisr2{2}(i)])/max(mean([thisr2{1}; thisr2{2}])), 'markersize', 10);
                hold on;
        end
%     plot(thisr2{2}, thisr2{1}, 'k.', 'markersize', 10); hold on
%     lsline;
    x = 0:.005:1; y=x;
    plot(x,y,'-', 'Color', [0.5,0.5,0.5]);
    
    x = thisr2{2};
%     x = [x(:) ones(size(x(:)))];
    y = thisr2{1};
%     beta = x\y; %((x'*x)^-1)*(x')*y'; beta = ((x'*w*x)^-1)*(x'*w)*y';
%     if ~isempty(beta)
%         plot([0 1], [0 1]*beta(1)+beta(2), 'r-');
%     end
    if length(x)>2
        mdl =fitlm(x,y);
        plot([0 1], [0 1]*double(mdl.Coefficients(2,1))+double(mdl.Coefficients(1,1)), '-', 'Color', [0.9 0.3 0.6]);
        title(sprintf('%s r^2\nslope:%0.2f (%s)\noffset:%0.2f (%s) <r2=%0.2f>', [hemi{h} ROIs{roi}], double(mdl.Coefficients(2,1)), pvaldisp(double(mdl.Coefficients(2,4))), ...
        double(mdl.Coefficients(1,1)), pvaldisp(double(mdl.Coefficients(1,4))), mdl.Rsquared.Ordinary));
    else
        title(sprintf('%s r^2',[hemi{h} ROIs{roi}]));
    end
    maxr2 = max([thisr2{2} thisr2{1}]);
    xaxis(r2cutoff, round(maxr2*10)/10);
    yaxis(r2cutoff, round(maxr2*10)/10);
    axis square;
    xlabel('r^2: BarsTaskFixation');
    ylabel('r^2: BarsTask');
    
weights = mean([thisr2{1}; thisr2{2}]);

 % plot RF center (x)
    subplot(roiGroupN,ncols,2+(roi-subtract)*ncols)
    for i = 1:length(thisX{2})
%       if ~isnan(thisr2{2}(i)) && ~isnan(thisr2{1}(i))
        plot(thisX{2}(i),thisX{1}(i),'k.','Color',1-c*mean([thisr2{1}(i) thisr2{2}(i)])/max(mean([thisr2{1}; thisr2{2}])), 'markersize', 10);
        hold on;
%       end
    end
%     plot(thisX{2}, thisX{1}, 'k.', 'markersize', 10); hold on
%     lsline;
    x = minX(2):.5:maxX(2); y=x;
    plot(x,y,'-', 'Color', [0.5,0.5,0.5]);
        x = thisX{2};
        y = thisX{1};

    if length(x)>2
        mdl =fitlm(x,y);
        plot([minX(1) maxX(1)], [minX(1) maxX(1)]*double(mdl.Coefficients(2,1))+double(mdl.Coefficients(1,1)), ':', 'Color', [0.9 0.6 0.7]);

        mdl =fitlm(x,y,'linear','Weights',weights);
        plot([minX(1) maxX(1)], [minX(1) maxX(1)]*double(mdl.Coefficients(2,1))+double(mdl.Coefficients(1,1)), '-', 'Color', [0.9 0.3 0.6]);
        title(sprintf('%s RF center (x)\nslope:%0.2f (%s)\noffset:%0.2f (%s) <r2=%0.2f>', [hemi{h} ROIs{roi}], double(mdl.Coefficients(2,1)), pvaldisp(double(mdl.Coefficients(2,4))), ...
        double(mdl.Coefficients(1,1)), pvaldisp(double(mdl.Coefficients(1,4))), mdl.Rsquared.Ordinary));
    else
        title(sprintf('%s RF center (x)',[hemi{h} ROIs{roi}]));
    end
    
    xaxis(min(minX(2),minX(1)), max(maxX(1),maxX(2)));
    yaxis(min(minX(1),minX(2)), max(maxX(1),maxX(2)));
    axis square;
    xlabel('x: BarsTaskFixation');
    ylabel('x: BarsTask');
%     title(sprintf('%s RF center (x)',[hemi{h} ROIs{roi}]));
    
 % plot RF center (y)
    subplot(roiGroupN,ncols,3+(roi-subtract)*ncols);
    for i = 1:length(thisY{2})
%       if ~isnan(thisr2{2}(i)) && ~isnan(thisr2{1}(i))
        plot(thisY{2}(i),thisY{1}(i),'k.','Color',1-c*mean([thisr2{1}(i) thisr2{2}(i)])/max(mean([thisr2{1}; thisr2{2}])), 'markersize', 10);
        hold on;
%       end
    end
%     plot(thisY{2}, thisY{1}, 'k.', 'markersize', 10); hold on
%     lsline;
    x = minY(2):.5:maxY(2); y=x;
    plot(x,y,'-', 'Color', [0.5,0.5,0.5]);
    
    x = thisY{2};
        y = thisY{1};

    if length(x)>2
        mdl =fitlm(x,y);
        plot([minY(1) maxY(1)], [minY(1) maxY(1)]*double(mdl.Coefficients(2,1))+double(mdl.Coefficients(1,1)), ':', 'Color', [0.9 0.6 0.7]);

        mdl =fitlm(x,y,'linear','Weights',weights);
        plot([minY(1) maxY(1)], [minY(1) maxY(1)]*double(mdl.Coefficients(2,1))+double(mdl.Coefficients(1,1)), '-', 'Color', [0.9 0.3 0.6]);
        title(sprintf('%s RF center (y)\nslope:%0.2f (%s)\noffset:%0.2f (%s) <r2=%0.2f>', [hemi{h} ROIs{roi}], double(mdl.Coefficients(2,1)), pvaldisp(double(mdl.Coefficients(2,4))), ...
        double(mdl.Coefficients(1,1)), pvaldisp(double(mdl.Coefficients(1,4))), mdl.Rsquared.Ordinary));
    else
        title(sprintf('%s RF center (y)',[hemi{h} ROIs{roi}]));
    end
    xaxis(min(minY(2),minY(1)), max(maxY(1),maxY(2)));
    yaxis(min(minY(1),minY(2)), max(maxY(1),maxY(2)));
    axis square;
    xlabel('y: BarsTaskFixation');
    ylabel('y: BarsTask');
%     title(sprintf('%s RF center (y)',[hemi{h} ROIs{roi}]));
    
 % plot rfHalfWidth
    subplot(roiGroupN,ncols,4+(roi-subtract)*ncols)
    for i = 1:length(thisRfHalfWidth{2})
%       if ~isnan(thisr2{2}(i)) && ~isnan(thisr2{1}(i))
        plot(thisRfHalfWidth{2}(i),thisRfHalfWidth{1}(i),'k.','Color',1-c*mean([thisr2{1}(i) thisr2{2}(i)])/max(mean([thisr2{1}; thisr2{2}])), 'markersize', 10);
        hold on;
%       end
    end
%     plot(thisRfHalfWidth{2}, thisRfHalfWidth{1}, 'k.', 'markersize', 10); hold on
%     lsline
    x = min([thisRfHalfWidth{2} thisRfHalfWidth{1}]):.5:max([thisRfHalfWidth{2} thisRfHalfWidth{1}]); y=x;
    plot(x,y,'-', 'Color', [0.5,0.5,0.5]);
    
    x = thisRfHalfWidth{2};
    y = thisRfHalfWidth{1};

    if length(x)>2
         mdl =fitlm(x,y);
        plot([0 max([thisRfHalfWidth{1} thisRfHalfWidth{2}])], [0 max([thisRfHalfWidth{1} thisRfHalfWidth{2}])]*double(mdl.Coefficients(2,1))+double(mdl.Coefficients(1,1)), ':', 'Color', [0.9 0.6 0.7]);

        mdl =fitlm(x,y,'linear','Weights',weights);
        plot([0 max([thisRfHalfWidth{1} thisRfHalfWidth{2}])], [0 max([thisRfHalfWidth{1} thisRfHalfWidth{2}])]*double(mdl.Coefficients(2,1))+double(mdl.Coefficients(1,1)), '-', 'Color', [0.9 0.3 0.6]);
        title(sprintf('%s RF half width\nslope:%0.2f (%s)\noffset:%0.2f (%s) <r2=%0.2f>', [hemi{h} ROIs{roi}], double(mdl.Coefficients(2,1)), pvaldisp(double(mdl.Coefficients(2,4))), ...
        double(mdl.Coefficients(1,1)), pvaldisp(double(mdl.Coefficients(1,4))), mdl.Rsquared.Ordinary));
    else
        title(sprintf('%s RF half width', [hemi{h} ROIs{roi}]));
    end
    
     xaxis(min([thisRfHalfWidth{1} thisRfHalfWidth{2}]), max([thisRfHalfWidth{1} thisRfHalfWidth{2}]));
     yaxis(min([thisRfHalfWidth{1} thisRfHalfWidth{2}]), max([thisRfHalfWidth{1} thisRfHalfWidth{2}]));
%      xaxis(0,20)
%      yaxis(0,20)
    axis square;
    xlabel('rfHalfWidth: BarsTaskFixation');
    ylabel('rfHalfWidth: BarsTask');
    
    
 % plot eccentricity
    subplot(roiGroupN,ncols,5+(roi-subtract)*ncols)
    for i = 1:length(thisEccentricity{2})
%       if ~isnan(thisr2{2}(i)) && ~isnan(thisr2{1}(i))
        plot(thisEccentricity{2}(i),thisEccentricity{1}(i),'k.','Color',1-c*mean([thisr2{1}(i) thisr2{2}(i)])/max(mean([thisr2{1}; thisr2{2}])), 'markersize', 10);
        hold on;
%       end
    end
%     plot(thisEccentricity{2}, thisEccentricity{1}, 'k.', 'markersize', 10); hold on
%     lsline
    x = min([thisEccentricity{2} thisEccentricity{1}]):1:max([thisEccentricity{2} thisEccentricity{1}]); y=x;
    plot(x,y,'-', 'Color', [0.5,0.5,0.5]);
    x = thisEccentricity{2};
    y = thisEccentricity{1};

    if length(x)>2
        mdl = fitlm(x,y);
        plot([0 max([thisEccentricity{1} thisEccentricity{2}])], [0 max([thisEccentricity{1} thisEccentricity{2}])]*double(mdl.Coefficients(2,1))+double(mdl.Coefficients(1,1)), ':', 'Color', [0.9 0.6 0.7]);
        mdl =fitlm(x,y,'linear','Weights',weights);
        plot([0 max([thisEccentricity{1} thisEccentricity{2}])], [0 max([thisEccentricity{1} thisEccentricity{2}])]*double(mdl.Coefficients(2,1))+double(mdl.Coefficients(1,1)), '-', 'Color', [0.9 0.3 0.6]);
        title(sprintf('%s Eccentricity\nslope:%0.2f (%s)\noffset:%0.2f (%s) <r2=%0.2f>', [hemi{h} ROIs{roi}], double(mdl.Coefficients(2,1)), pvaldisp(double(mdl.Coefficients(2,4))), ...
        double(mdl.Coefficients(1,1)), pvaldisp(double(mdl.Coefficients(1,4))), mdl.Rsquared.Ordinary));
    else
        title(sprintf('%s Eccentricity',[hemi{h} ROIs{roi}]));
    end
    
    xaxis(min([thisEccentricity{1} thisEccentricity{2}]), max([thisEccentricity{1} thisEccentricity{2}]));
    yaxis(min([thisEccentricity{1} thisEccentricity{2}]), max([thisEccentricity{1} thisEccentricity{2}]));
%      xaxis(0,50)
%      yaxis(0,50)
    axis square;
    xlabel('Eccentricity: BarsTaskFixation');
    ylabel('Eccentricity: BarsTask');
end


