function [r2, x, y, Eccentricity, Size] = createMat4R
load('bothhemi_IPS5.mat')

roiList = {'V1','V2v','V2d','V3v','V3d','V4','V3A','V3B','LO1','LO2','IPS0','IPS1','IPS2','IPS3','IPS4','IPS5'};
roiName = {'V1','V2','V3','V4','V3A','V3B','LO1','IPS0','IPS1','IPS2','IPS3','IPS4','IPS5'};

for roi = 1:length(roiName)

    if ~(strcmp(roiName{roi},'V2') || strcmp(roiName{roi},'V3'))
        
        roiInd = find(strcmp(roiList, roiName{roi}));

        for cond = 1:2
            %NUMBER IN ROI NAME!!!!!!! NOT the whole roiList
            r2(roi*2+cond-2,1) = roi; %ROI number
            r2(roi*2+cond-2,2) = cond; %condition num (1 for task, 2 for fixation)
            r2(roi*2+cond-2,3) = mean(bothHemi.thisr2{roiInd}{cond}); %mean for each cond
            r2(roi*2+cond-2,4) = std(bothHemi.thisr2{roiInd}{cond})/sqrt(length(bothHemi.thisr2{roiInd}{cond})); % standard error
            
            x(roi*2+cond-2,1) = roi; %ROI number
            x(roi*2+cond-2,2) = cond; %condition num (1 for task, 2 for fixation)
            x(roi*2+cond-2,3) = mean(bothHemi.thisX{roiInd}{cond}); %mean for each cond
            x(roi*2+cond-2,4) = std(bothHemi.thisX{roiInd}{cond})/sqrt(length(bothHemi.thisX{roiInd}{cond})); % standard error
            
            y(roi*2+cond-2,1) = roi; %ROI number
            y(roi*2+cond-2,2) = cond; %condition num (1 for task, 2 for fixation)
            y(roi*2+cond-2,3) = mean(bothHemi.thisY{roiInd}{cond}); %mean for each cond
            y(roi*2+cond-2,4) = std(bothHemi.thisY{roiInd}{cond})/sqrt(length(bothHemi.thisY{roiInd}{cond})); % standard error

            Eccentricity(roi*2+cond-2,1) = roi; %ROI number
            Eccentricity(roi*2+cond-2,2) = cond; %condition num (1 for task, 2 for fixation)
            Eccentricity(roi*2+cond-2,3) = mean(bothHemi.thisEccentricity{roiInd}{cond}); %mean for each cond
            Eccentricity(roi*2+cond-2,4) = std(bothHemi.thisEccentricity{roiInd}{cond})/sqrt(length(bothHemi.thisEccentricity{roiInd}{cond})); % standard error
            
            Size(roi*2+cond-2,1) = roi; %ROI number
            Size(roi*2+cond-2,2) = cond; %condition num (1 for task, 2 for fixation)
            Size(roi*2+cond-2,3) = mean(bothHemi.thisRfHalfWidth{roiInd}{cond}); %mean for each cond
            Size(roi*2+cond-2,4) = std(bothHemi.thisRfHalfWidth{roiInd}{cond})/sqrt(length(bothHemi.thisRfHalfWidth{roiInd}{cond})); % standard error
        end
        
    else %V2, V3
        if strcmp(roiName{roi},'V2')
            ind1 = 2; ind2 = 3;
        elseif strcmp(roiName{roi},'V3')
            ind1 = 4; ind2 = 5;
        end
            
        for cond = 1:2
            thisr2{cond} = [bothHemi.thisr2{ind1}{cond}, bothHemi.thisr2{ind2}{cond}];
            thisX{cond} = [bothHemi.thisX{ind1}{cond}, bothHemi.thisX{ind2}{cond}];
            thisY{cond} = [bothHemi.thisY{ind1}{cond}, bothHemi.thisY{ind2}{cond}];
            thisEccentricity{cond} = [bothHemi.thisEccentricity{ind1}{cond}, bothHemi.thisEccentricity{ind2}{cond}];
            thisRfHalfWidth{cond} = [bothHemi.thisRfHalfWidth{ind1}{cond}, bothHemi.thisRfHalfWidth{ind2}{cond}];
        end
        for cond= 1:2
            r2(roi*2+cond-2,1) = roi; %ROI number
            r2(roi*2+cond-2,2) = cond; %condition num (1 for task, 2 for fixation)
            r2(roi*2+cond-2,3) = mean(thisr2{cond}); %mean for each cond
            r2(roi*2+cond-2,4) = std(thisr2{cond})/sqrt(length(thisr2{cond})); % standard error
            
            x(roi*2+cond-2,1) = roi; %ROI number
            x(roi*2+cond-2,2) = cond; %condition num (1 for task, 2 for fixation)
            x(roi*2+cond-2,3) = mean(thisX{cond}); %mean for each cond
            x(roi*2+cond-2,4) = std(thisX{cond})/sqrt(length(thisX{cond})); % standard error
            
            y(roi*2+cond-2,1) = roi; %ROI number
            y(roi*2+cond-2,2) = cond; %condition num (1 for task, 2 for fixation)
            y(roi*2+cond-2,3) = mean(thisY{cond}); %mean for each cond
            y(roi*2+cond-2,4) = std(thisY{cond})/sqrt(length(thisY{cond})); % standard error

            Eccentricity(roi*2+cond-2,1) = roi; %ROI number
            Eccentricity(roi*2+cond-2,2) = cond; %condition num (1 for task, 2 for fixation)
            Eccentricity(roi*2+cond-2,3) = mean(thisEccentricity{cond}); %mean for each cond
            Eccentricity(roi*2+cond-2,4) = std(thisEccentricity{cond})/sqrt(length(thisEccentricity{cond})); % standard error
            
            Size(roi*2+cond-2,1) = roi; %ROI number
            Size(roi*2+cond-2,2) = cond; %condition num (1 for task, 2 for fixation)
            Size(roi*2+cond-2,3) = mean(thisRfHalfWidth{cond}); %mean for each cond
            Size(roi*2+cond-2,4) = std(thisRfHalfWidth{cond})/sqrt(length(thisRfHalfWidth{cond})); % standard error
        end
    end
end


      
