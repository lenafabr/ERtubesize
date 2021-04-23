% find min distance between two curves

tubepath = CL.ROIgroups(1).tubeROIs(2).interppath;
sheetbound = CL.ROIgroups(1).sheetROI.bound;

mindist = inf;
for pc = 1:size(sheetbound,1)
    dists = sum((tubepath - sheetbound(pc,:)).^2,2);
    mindist = min(mindist,min(dists));
end


%% load in WT data
load('../results/400ngCtrl_b_FITC_Calnexin_membrane.mat')
erodeum=0.2; dilum =0.2;
RestimatesWTLaura = [];
for cc = 1:length(allcells)    
    CL = allcells(cc);
    CL.reprocessROIs(erodeum,dilum);
    for sc = 1:length(CL.ROIgroups)       
        CL.ROIgroups(sc).Restimate = getRestimates(CL,CL.ROIgroups(sc),struct('mintubelen',5));
        RestimatesWTLaura(end+1) = CL.ROIgroups(sc).Restimate;
    end    
end

%% more WT data
load('../results/200ngCtrl_Calnexin_membrane.mat')
erodeum=0.2; dilum =0.2;
Restimates200WTLaura = [];
for cc = 1:length(allcells)    
    CL = allcells(cc);
    CL.reprocessROIs(erodeum,dilum);
    for sc = 1:length(CL.ROIgroups)       
        CL.ROIgroups(sc).Restimate = getRestimates(CL,CL.ROIgroups(sc),struct('mintubelen',5));
        Restimates200WTLaura(end+1) = CL.ROIgroups(sc).Restimate;
    end    
end

%% RTN4OE data
load('../results/200ngRTN4OE_Calnexin_membrane.mat')
erodeum=0.2; dilum =0.2;
Restimates200RTN4OE = [];
for cc = 1:length(allcells)    
    CL = allcells(cc);
    CL.reprocessROIs(erodeum,dilum);
    for sc = 1:length(CL.ROIgroups)       
        CL.ROIgroups(sc).Restimate = getRestimates(CL,CL.ROIgroups(sc),struct('mintubelen',5));
        Restimates200RTN4OE(end+1) = CL.ROIgroups(sc).Restimate;
    end    
end