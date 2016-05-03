function v = combineROIsTaskFix()

v = getMLRView();
roiList = {'V1';'V2d';'V3d';'V2v';'V3v';'V4';'V3A';'V3B';'LO1';'LO2';'MT';'IPS0';'IPS1';'IPS2';'IPS3';'IPS4';'IPS5'; 'SPL1'};


hemi = {'l', 'r'};

disp(sprintf('(blabla) Running --- on'))
disp(sprintf('%s\t', roiList{:}))
allROIs = askuser('Do you wish to continue?',1);
if ~allROIs
    return;
end

for h = 1:length(hemi)

    for r = 1:length(roiList)
        
        roi1 = viewGet(v,'roiNum',[hemi{h} roiList{r} '_t']);
        roi2 = viewGet(v,'roiNum', [hemi{h} roiList{r} '_f']);
        if isempty(roi1) || isempty(roi2)
            continue
        end
        v = combineROIs(v,roi1,roi2,'Intersection', [hemi{h} roiList{r} '_C']);
    end
end
        
        