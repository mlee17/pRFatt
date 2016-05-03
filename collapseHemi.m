function [hemiCollapsed, df] = collapseHemi(ROIStruct)

if ieNotDefined('ROIStruct')
%     ROIStruct = '/data/mglretinotopy/s316RoiStruct_c.mat';
    return
end

sbj = ROIStruct{1}{1}.Subject;

% hemiCollapsed = struct('Task', struct( cell(1, length(ROIStruct{1})), 'Fix', cell(1, length(ROIStruct{1})));
% hemiCollapsed.
tasktype = {'Task','Fix'};
for roi = 1:length(ROIStruct{1})
    
    for hemi = 1:2
        if isempty(ROIStruct{hemi}{roi})
            for type = 1:2
                ROIStruct{hemi}{roi}(type).roiName = [];
                ROIStruct{hemi}{roi}(type).thisr2 = []; ROIStruct{hemi}{roi}(type).thisX = []; ROIStruct{hemi}{roi}(type).thisY = []; 
                ROIStruct{hemi}{roi}(type).thisEccentricity = []; ROIStruct{hemi}{roi}(type).thisRfHalfWidth=[]; 
            end
        end
        
        
    for type = 1:2
        hemiCollapsed.(tasktype{type}){roi}.r2 = [ROIStruct{1}{roi}(type).thisr2 ROIStruct{2}{roi}(type).thisr2];
        hemiCollapsed.(tasktype{type}){roi}.x = [ROIStruct{1}{roi}(type).thisX ROIStruct{2}{roi}(type).thisX];
        hemiCollapsed.(tasktype{type}){roi}.y = [ROIStruct{1}{roi}(type).thisY ROIStruct{2}{roi}(type).thisY];
        hemiCollapsed.(tasktype{type}){roi}.eccentricity = [ROIStruct{1}{roi}(type).thisEccentricity ROIStruct{2}{roi}(type).thisEccentricity];
        hemiCollapsed.(tasktype{type}){roi}.width = [ROIStruct{1}{roi}(type).thisRfHalfWidth ROIStruct{2}{roi}(type).thisRfHalfWidth];
    end
    end
    
    if ~isempty(ROIStruct{1}{roi}(1).roiName)
        roiList{roi} = ROIStruct{1}{roi}(1).roiName(2:end);
    elseif ~isempty(ROIStruct{2}{roi}(1).roiName)
         roiList{roi} = ROIStruct{2}{roi}(1).roiName(2:end);
    else
        roiList{roi}=[];
    end

end

df = data2Rformat(hemiCollapsed,roiList);

%     df_r2 = data2Rformat(hemiCollapsed,roiList, 'r2');
%     df_x = data2Rformat(hemiCollapsed,roiList, 'x');
%     df_y = data2Rformat(hemiCollapsed,roiList, 'y');
%     df_eccentricity = data2Rformat(hemiCollapsed,roiList, 'eccentricity');
%     df_width = data2Rformat(hemiCollapsed,roiList, 'width');
          

%%
function df = data2Rformat(structure, roiList)
varList = {'r2','x','y','eccentricity','width'};
for r = 1:length(roiList)
    for whichparam = 1:length(varList)
        df.(roiList{r}(1:end-2)).(varList{whichparam}) = [structure.Task{r}.(varList{whichparam})'; structure.Fix{r}.(varList{whichparam})'];
    end
    
    len = length(structure.Task{r}.(varList{whichparam}));
    df.(roiList{r}(1:end-2)).cond(1:len,1) = 1;
    df.(roiList{r}(1:end-2)).cond(len+1:2*len,1) = 2;
    
end


return
%         df.(roiList{roi}) = struct(sprintf(varList{whichparam}), [structure.Task{roi}.(varList{whichparam})'; structure.Fix{roi}.(varList{whichparam})'])
            
    

% 
%       df.(roiList{roi})(:).(varList{whichparam}) = [structure.Task{roi}.(varList{whichparam})'; structure.Fix{roi}.(varList{whichparam})'];
%         len = length(structure.Task{roi}.(varList{whichparam}));
%         df.(roiList{roi})(1:len,2).cond = 1;
%         
%         df.(roiList{roi})(len+1,2) = 2;
%     
%     
%     df.(roiList{roi}) = [hemiCollapsed.Task{roi}.(whichparam)'; hemiCollapsed.Fix{roi}.(whichparam)'];
%     df.(roiList{roi})(1:length(df.(roiList{roi}))/2,2) = 1;
%     df.(roiList{roi})((length(df.(roiList{roi}))/2)+1,2) = 2;


% 
% function df = data2Rformat(structure, roiList, whichparam)
% for roi = 1:length(roiList)
%     df.(roiList{roi}(1:end-2)) = [structure.Task{roi}.(whichparam)'; structure.Fix{roi}.(whichparam)'];
%     df.(roiList{roi}(1:end-2))(1:length(df.(roiList{roi}))/2,2) = 1;
%     df.(roiList{roi}(1:end-2))((length(df.(roiList{roi}))/2)+1,2) = 2;
% end
% 
% 

 
    
