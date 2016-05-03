
%roiname = {'V1','V2v','V2d','V3v','V3d','V3A','V3B','LO1','LO2','V4','IPS0','IPS1','IPS2','IPS3','IPS4'};
roiList = {'V1','V2v','V2d','V3v','V3d','V4','V3A','V3B','LO1','LO2','IPS0','IPS1','IPS2','IPS3','IPS4'};
roiName = {'V1','V2','V3','V4','V3A','V3B','LO1','IPS0','IPS1','IPS2','IPS3','IPS4'};

varList = {'r2','x','y','eccentricity','width'};

for r = 1:length(roiName)
    for param = 1:length(varList)
        if ~(strcmp(roiName{r},'V2') || strcmp(roiName{r},'V3'))
            for c = 1:2
                barPrep.(varList{param})(r,c) = mean(df.(roiName{r}).(varList{param})(df.(roiName{r}).cond ==c));
%             barPrep.(varList{parmam})(r,2) = mean(df.(roiName{r}).(varList{parmam})(df.(roiName{r}).cond ==2));
        
                errors.(varList{param})(r,c) = std(df.(roiName{r}).(varList{param})(df.(roiName{r}).cond ==c))/sqrt(sum(df.(roiName{r}).cond ==1));
%             errors.(varList{parmam})(r,2) = std(df.(roiName{r}).(varList{parmam})(df.(roiName{r}).cond ==2))/sqrt(sum(df.(roiName{r}).cond ==1));
            end
        else
            for c=1:2
                tmpData = [df.([roiName{r} 'd']).(varList{param})(df.([roiName{r} 'd']).cond ==c); df.([roiName{r} 'v']).(varList{param})(df.([roiName{r} 'v']).cond ==c)];
                barPrep.(varList{param})(r,c) = mean(tmpData);
                errors.(varList{param})(r,c) = std(tmpData) / sqrt(length(tmpData));
            end
        end
    end
end


% Usage: handles = barweb(barvalues, errors, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, error_sides, legend_type)
hd1 = barweb(barPrep.width, errors.width,[], roiName, [],[],'pRF size (degrees)', winter, [], {'task','fixation'}, 2, 'plot');

hd2 = barweb(barPrep.width, errors.width,[], roiName, [],[],'pRF size (degrees)', winter, [], {'task','fixation'}, 2, 'plot');
            
