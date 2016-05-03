% numFit.m
%
%        $Id:$ 
%      usage: numFit()
%         by: Minyoung Lee / Justin Gardner
%       date: 12/15/2015
%    purpose: Fit pRF model to voxel
%
function retval = pRFMinyoung(v,overlayNum,scanNum,x,y,s,roi)

% get time series
tSeries = squeeze(loadTSeries(v,scanNum,s,[],x,y));tSeries = tSeries(:);

% percent signal change
tSeries = (tSeries-mean(tSeries))/mean(tSeries);

% get basic paramters (like for canonical HDR)
[v fitParams] = pRF(v,[],'justGetParams=1','defaultParams=1','scanList',scanNum);
fitParams = fitParams.pRFFit;

% set the inital parameters to start off the model with
fitParams.initParams = [5.63 -4.78 3.09];

% set algorithm for nonlinear optimization and some parameters
fitParams.algorithm = 'nelder-mead';
fitParams.optimParams = optimset('MaxIter',inf);
fitParams.framePeriod = viewGet(v,'framePeriod',scanNum);

% set concatInfo
fitParams.concatInfo = getConcatInfo(v,scanNum);

% get this to the stimulus sequence
fitParams.stim = getStim(v,scanNum,fitParams);

% display the stimulus (just for debugging)
%mlrVol(fitParams.stim{1}.im)

% Just for testing
r = getModelResidual(fitParams.initParams,tSeries,fitParams);

% find best parameters
[params fval exitflag] = fminsearch(@getModelResidual,fitParams.initParams,fitParams.optimParams,tSeries,fitParams);
r

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getModelResidual   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [residual modelResponse rfModel r amp] = getModelResidual(params,tSeries,fitParams,justGetModel)

residual = [];
if nargin < 4, justGetModel = 0;end

% get the model response
% convert parameter array into a parameter strucutre
p = getFitParams(params,fitParams);

% FIX, compute the gaussian receptive field given the fitParams
rfModel = getRFModel(p,fitParams);
%rfModel = exp(-(((fitParams.stim{1}.x - p.x).^2 + (fitParams.stim{1}.y - p.y).^2) / (2*(p.std)^2)));

% FIX, compute the RF response to the stimulus (using fitParams.stim which should have the stimulus sequence - stim.im are the images. The other fields contain the units for stim.im: stim.x, stim.y contain what the x an y position in degrees is of each point in stim.im and stim.t contains the time in seconds)
% modelResponse = zeros(1,length(tSeries));
%modelResponse = computeModelResponse(rfModel, fitParams.stim,...)
% for t = 1:length(tSeries)
%     modelResponse(t) = sum(sum(fitParams.stim{1}.im(:,:,t) .* rfModel));
% end
modelResponse = computeModelResponse(rfModel, fitParams.stim{1});

% get a model hrf
hrf = getCanonicalHRF(p.canonical,fitParams.framePeriod);

% and convolve in time.
modelResponse = convolveModelResponseWithHRF(modelResponse,hrf);

% with no filtering, just remove mean
modelResponse = modelResponse - mean(modelResponse);

% scale the model response
amp = pinv(modelResponse') * tSeries;
modelResponse = modelResponse * amp;

% compute correlation of model with time series
r = corr(tSeries(:),modelResponse(:));


dispFit = 1;
if dispFit
  mlrSmartfig('numFit','reuse');clf;
  subplot(4,4,[1:3 5:7 9:11 13:15]);
  plot(100*tSeries,'k.-');hold on;
  plot(100*modelResponse,'r-');
  xlabel('Time (vols)');
  ylabel('BOLD (% signal change)');
  title(sprintf('x: %f y: %f width: %f',p.x,p.y,p.std));
  subplot(4,4,8);
  plot(hrf.time, hrf.hrf);
  subplot(4,4,12);
  imagesc(fitParams.stim{1}.x(:,1), fitParams.stim{1}.y(1,:), flipud(rfModel'))
  drawnow
end


% for nelder-mead just compute correlation and return negative
% since nelder-mean is optimizing for the least value
if strcmp(lower(fitParams.algorithm),'nelder-mead')
  residual = -r;
end
%%%%%%%%%%%%%%%%%%%%%%
%%    getRFModel    %%
function rfModel = getRFModel(p, fitParams)
rfModel = exp(-(((fitParams.stim{1}.x - p.x).^2)/(2*(p.std^2))+((fitParams.stim{1}.y - p.y).^2)/(2*(p.std^2))));
%%%%%%%%%%%%%%%%%%%%%%
%%   computeModelResponse   %%
function modelResponse = computeModelResponse(rfModel, stim)
modelResponse = zeros(1,length(stim.t));
for t = 1:length(stim.t)
    modelResponse(t) = sum(sum(stim.im(:,:,t) .* rfModel));
end
%%%%%%%%%%%%%%%%%%%%%%
%%   getFitParams   %%
%%%%%%%%%%%%%%%%%%%%%%
function p = getFitParams(params,fitParams)

p.x = params(1);
p.y = params(2);
p.std = params(3);

% use a fixed single gaussian
p.canonical.type = 'gamma';
p.canonical.lengthInSeconds = 25;
p.canonical.timelag = fitParams.timelag;
p.canonical.tau = fitParams.tau;
p.canonical.exponent = fitParams.exponent;
p.canonical.offset = 0;
p.canonical.diffOfGamma = fitParams.diffOfGamma;
p.canonical.amplitudeRatio = fitParams.amplitudeRatio;
p.canonical.timelag2 = fitParams.timelag2;
p.canonical.tau2 = fitParams.tau2;
p.canonical.exponent2 = fitParams.exponent2;
p.canonical.offset2 = 0;

%%%%%%%%%%%%%%%%%%%%%
%%   getGammaHRF   %%
%%%%%%%%%%%%%%%%%%%%%
function fun = getGammaHRF(time,p)

fun = thisGamma(time,1,p.timelag,p.offset,p.tau,p.exponent)/100;
% add second gamma if this is a difference of gammas fit
if p.diffOfGamma
  fun = fun - thisGamma(time,p.amplitudeRatio,p.timelag2,p.offset2,p.tau2,p.exponent2)/100;
end

%%%%%%%%%%%%%%%%%%%
%%   thisGamma   %%
%%%%%%%%%%%%%%%%%%%
function gammafun = thisGamma(time,amplitude,timelag,offset,tau,exponent)

exponent = round(exponent);
% gamma function
gammafun = (((time-timelag)/tau).^(exponent-1).*exp(-(time-timelag)/tau))./(tau*factorial(exponent-1));

% negative values of time are set to zero,
% so that the function always starts at zero
gammafun(find((time-timelag) < 0)) = 0;

% normalize the amplitude
if (max(gammafun)-min(gammafun))~=0
  gammafun = (gammafun-min(gammafun)) ./ (max(gammafun)-min(gammafun));
end
gammafun = (amplitude*gammafun+offset);


%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getCanonicalHRF   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function hrf = getCanonicalHRF(params,sampleRate)

hrf.time = 0:sampleRate:params.lengthInSeconds;
hrf.hrf = getGammaHRF(hrf.time,params);

% normalize to amplitude of 1
hrf.hrf = hrf.hrf / max(hrf.hrf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   convolveModelResponseWithHRF   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function modelTimecourse = convolveModelResponseWithHRF(modelTimecourse,hrf)

n = length(modelTimecourse);
modelTimecourse = conv(modelTimecourse,hrf.hrf);
modelTimecourse = modelTimecourse(1:n);

%%%%%%%%%%%%%%%%%
%    getStim    %
%%%%%%%%%%%%%%%%%
function stim = getStim(v,scanNum,fitParams)

% get stimfile
stimfile = viewGet(v,'stimfile',scanNum);
% get volume to trigger ratio
volTrigRatio = viewGet(v,'auxParam','volTrigRatio',scanNum);
% check if global matches
groupNum = viewGet(v,'curGroup');
global gpRFFitStimImage
if (isfield(fitParams,'recomputeStimImage') && fitParams.recomputeStimImage) || isempty(gpRFFitStimImage) || (gpRFFitStimImage.scanNum ~= scanNum)  || (gpRFFitStimImage.groupNum ~= groupNum) || (gpRFFitStimImage.xFlip ~= fitParams.xFlipStimulus) || (gpRFFitStimImage.yFlip ~= fitParams.yFlipStimulus) || (gpRFFitStimImage.timeShift ~= fitParams.timeShiftStimulus)
  disp(sprintf('(pRFFit) Computing stim image'));
  % if no save stim then create one
  stim = pRFGetStimImageFromStimfile(stimfile,'volTrigRatio',volTrigRatio,'xFlip',fitParams.xFlipStimulus,'yFlip',fitParams.yFlipStimulus,'timeShift',fitParams.timeShiftStimulus,'verbose',fitParams.verbose,'saveStimImage',fitParams.saveStimImage,'recomputeStimImage',fitParams.recomputeStimImage);
  % check for averages
  stim = checkStimForAverages(v,scanNum,viewGet(v,'curGroup'),stim,fitParams.concatInfo,fitParams.stimImageDiffTolerance);
  if isempty(stim),return,end
  % make into cell array
  stim = cellArray(stim);
  % save stim image in global
  gpRFFitStimImage.scanNum = scanNum;
  gpRFFitStimImage.groupNum = groupNum;
  gpRFFitStimImage.xFlip = fitParams.xFlipStimulus;
  gpRFFitStimImage.yFlip = fitParams.yFlipStimulus;
  gpRFFitStimImage.timeShift = fitParams.timeShiftStimulus;
  gpRFFitStimImage.stim = stim;
else
  % otherwise load from global
  disp(sprintf('(pRFFit) Using precomputed stim image'));
  stim = gpRFFitStimImage.stim;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    checkStimForAverages    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stim ignoreMismatchStimfiles] = checkStimForAverages(v,scanNum,groupNum,stim,concatInfo,stimImageDiffTolerance)

ignoreMismatchStimfiles = false;  

% this function will check for some bad casses (like concat of concats etc)
% it will also check that all the component scans of an average have the
% same stim image and warn if they do not. It will then replace the stim cell
% array for the average with a single stim file, so that processing
% can continue as normal for pRFFit

% if not a cell, then ok, return
if ~iscell(stim),return,end

% first check for bad shiftList or refverseLIst
p = viewGet(v,'params',scanNum,groupNum);
if isfield(p,'params') && isfield(p.params,'shiftList') && any(p.params.shiftList~=0)
  disp(sprintf('(pRFFit) Component scan %s:%i has a shiftList that is non-zero (%s). pRFFit does not handle non-zero shifts in averages.',viewGet(v,'groupName',groupNum),scanNum,mlrnum2str(p.params.shiftList)));
  keyboard
end
if isfield(p,'params') && isfield(p.params,'reverseList') && any(p.params.reverseList~=0)
  disp(sprintf('(pRFFit) Component scan %s:%i has a reverseList that is non-zero (%s). pRFFit does not handle time-reversed time series in averages.',viewGet(v,'groupName',groupNum),scanNum,mlrnum2str(p.params.shiftList)));
  keyboard
end

% if is a cell, check to see if this is a concat or not
if ~isempty(concatInfo) && (concatInfo.isConcat)
  % this is a concat, so check each one of the elements
  [originalScanNum originalGroupNum] = viewGet(v,'originalScanNum',scanNum,groupNum);
  for i = 1:length(stim)
    % get concatInfo for original scan
    concatInfo = viewGet(v,'concatInfo',originalScanNum(i),originalGroupNum(i));
    if ~isempty(concatInfo)
      disp(sprintf('(pRFFit:checkStimForAverages) Detected concatenation of concatenations. pRFFit not implemented yet to handle this'));
      stim = [];
      keyboard
      return;
    end
    % check this next scan
    [stim{i} ignoreMismatchStimfiles] = checkStimForAverages(v,originalScanNum(i),originalGroupNum(i),stim{i},concatInfo,stimImageDiffTolerance);
    % if user has accepted all then set stimImageDiffTOlerance to infinity
    if isinf(ignoreMismatchStimfiles),stimImageDiffTolerance = inf;end
    if isempty(stim{i}),stim = [];return,end
  end
else
  % this for orignals
  [originalScanNum originalGroupNum] = viewGet(v,'originalScanNum',scanNum,groupNum);
  % if it is an original than check each element
  if ~isempty(originalScanNum)
    % check that this is not an average of a concat
    for i = 1:length(stim)
      % get concatInfo for original scan
      concatInfo = viewGet(v,'concatInfo',originalScanNum(i),originalGroupNum(i));
      if ~isempty(concatInfo)
	disp(sprintf('(pRFFit:checkStimForAverages) Detected average of a concatenations. pRFFit not implemented yet to handle this'));
	keyboard
	stim = [];
	return;
      end
      % see if it is an average of an average
      originalOfOriginalScanNum = viewGet(v,'originalScanNum',originalScanNum(i),originalGroupNum(i));
      if length(originalOfOriginalScanNum) > 1
	disp(sprintf('(pRFFit:checkStimForAverages) Detected average of an average. pRFFit not implemented yet to handle this'));
	keyboard
	stim = [];
	return;
      end
    end
    % ok, not an average of a concatenation/average so check all the stim files 
    % and warn if there are any inconsistencies
    for i = 1:length(stim)
      if ~isequalwithequalnans(stim{1}.im,stim{i}.im)    
	dispHeader
	disp(sprintf('(pRFFit:checkStimForAverages) !!! Average for %s:%i component scan %i does not match stimulus for other scans. If you wish to continue then this will use the stimfile associated with the first scan in the average !!!',viewGet(v,'groupName',groupNum),scanNum,originalScanNum(i)));
	% display which volumes are different
	diffVols = [];
	for iVol = 1:size(stim{1}.im,3)
	  if ~isequalwithequalnans(stim{1}.im(:,:,iVol),stim{i}.im(:,:,iVol))
	    diffVols(end+1) = iVol;
	  end
	end
	disp(sprintf('(pRFFit) Stimulus files are different at %i of %i vols (%0.1f%%): %s',length(diffVols),size(stim{1}.im,3),100*length(diffVols)/size(stim{1}.im,3),num2str(diffVols)));
	if 100*(length(diffVols)/size(stim{1}.im,3)) < stimImageDiffTolerance
	  disp(sprintf('(pRFFit) This could be for minor timing inconsistencies, so igorning. Set stimImageDiffTolerance lower if you want to stop the code when this happens'));
	else
	  % ask user if they want to continue (only if there is a difference of more than 10 vols	  
	  ignoreMismatchStimfiles = askuser('Do you wish to continue',1);
	  if ~ignoreMismatchStimfiles
	    stim = [];
	    return;
	  end
	end
	dispHeader
      end
    end
    % if we passed the above, this is an average of identical
    % scans, so just keep the first stim image since they are all the same
    stim = stim{1};
  end
end

%%%%%%%%%%%%%%%%%%%%%%%
%    getConcatInfo    %
%%%%%%%%%%%%%%%%%%%%%%%
function concatInfo = getConcatInfo(v,scanNum)

concatInfo = viewGet(v,'concatInfo',scanNum);

% if there is no concatInfo, then make one that will
% treat the scan as a single scan
if isempty(concatInfo)
  nFrames = viewGet(v,'nFrames',scanNum);
  concatInfo.isConcat = false;
  concatInfo.n = 1;
  concatInfo.whichScan = ones(1,nFrames);
  concatInfo.whichVolume = 1:nFrames;
  concatInfo.runTransition = [1 nFrames];
  concatInfo.totalJunkedFrames = viewGet(v,'totalJunkedFrames',scanNum);
  if length(concatInfo.totalJunkedFrames > 1)
    % first check for consistency in totalJunkedFrames
    if length(unique(concatInfo.totalJunkedFrames)) > 1
      disp(sprintf('(pRFFit) totalJunkedFrames are different for different members of component scans - could be an average in which different scans with different number of junked frames were removed. This could cause a problem in computing what the stimulus was for the average. The total junked frames count was: %s, but we will use %i as the actual value for computing the stimulus',num2str(fitParams.concatInfo.totalJunkedFrames),floor(median(fitParams.concatInfo.totalJunkedFrames))));
    end
    concatInfo.totalJunkedFrames = floor(median(concatInfo.totalJunkedFrames));
  end
else
  concatInfo.isConcat = true;
  if ~isfield(concatInfo,'totalJunkedFrames')
    concatInfo.totalJunkedFrames = viewGet(v,'totalJunkedFrames',scanNum);
  end
end
