
dirname = ['~/UCSD/data/Avezov/Tasuko/220421_COS7_Halo-Sec61b-TMR_for tubule width_live and fixed/'];
%
fglob = ['220421_COS7_Halo_Sec61b_TMR_fortubulewidth_live*.tif'];


files = dir([dirname fglob])

%% offset for dilation and erosion, to deal with resolution issues
resoffset = 0.2;

%%
cellct = 0;
%%
%cellct = 0;
for fc = 2%:length(files)
    fname = files(fc).name
    [filepath,name,ext] = fileparts(fname);
    
    CL = CellObjTubeSheet(name);
    memframe = 1; lumframe = 1;
    CL.loadCellData(dirname,fname,memframe,lumframe);        
    
    % show membrane and luminal images
    figure(1)    
    imshow(CL.imglum,[])
    title('luminal')
    
    %%
    ans = input('Use this cell? y/n \n','s')
    if (strcmp(ans,'y'))
        badcells(fc) = false;
    else
        badcells(fc) = true;
        continue
    end
    %% pick out background region
    figure(1)%; subplot(1,2,1)
    CL.getBackground(1);
            
    %%
    sheetcount = length(CL.ROIgroups)
    while 1        % sheets
        sheetcount = sheetcount + 1;
        
        %% pick out a sheet region
        figure(1)%; subplot(1,2,1)
        sheetROI = CL.getSheet(resoffset,1) 
        CL.ROIgroups(sheetcount).sheetROI = sheetROI;
        
        
        %% Pick out multiple tubules   
        tubect = length(CL.ROIgroups(sheetcount).tubeROIs);
        
        figure(1); %subplot(1,2,1)
        imshowpair(CL.imgmem,sheetROI.erodemask)
            
        while 1
            tubect = tubect+1;
            
            %% pick out a tube region
            
            tubeROI = CL.getTube(resoffset,1);
            
            if (tubect==1)
                CL.ROIgroups(sheetcount).tubeROIs = tubeROI;
            else
                CL.ROIgroups(sheetcount).tubeROIs(end+1) = tubeROI;
            end
            
            ans = input('More tubules? y/n \n','s')
            
            if (~strcmp(ans,'y'))
                break
            end
        end
        
        %% show masks
        CL.showROImasks(0)
        
                
        %% calculate radius estimate for this sheet and tubule group
        CL.ROIgroups(sheetcount).Restimate = CL.getRestimates(CL.ROIgroups(sheetcount));
        disp(sprintf('Group %d: Restimate %f', sheetcount,CL.ROIgroups(sheetcount).Restimate))
        
        %%
        ans = input('Pick another sheet? y/n \n','s');
        
        if (~strcmp(ans,'y'))
                break
        end
    end % done picking sheets
     
    %%
    cellct = cellct+1;
    allcells(cellct) = CL;
end

%%
save(['../results/220421_COS7_Halo_Sec61b_TMR_fortubulewidth_live.mat'])
%save('../results/210106_COS7_WT_Sec61_Halo_OG_ERmCherry.mat')
%save('210106_COS7_RTN4_KO_Sec61_Halo_OG_ER_mCherry.mat')
%save('../results/210223_COS7_WT_Sec61_Halo_OG_ER_mCherry.mat')
%save('../results/210223_COS7_RTN4_KO_Sec61_Halo_OG_ER_mCherry.mat')
%save('../results/210223_SHSY5Y_RTN4_KO_Sec61_Halo_OG_ER_mCherry.mat')
%% get averages
Restimates = [];
for cc = 1:length(allcells)    
    CL = allcells(cc);
    
    % exclude low resolution cells
    %lowres = contains(CL.Name,'_1');
    %if (lowres); continue; end
    
    for sc = 1:length(CL.ROIgroups)        
        Restimates(end+1) = CL.ROIgroups(sc).Restimate;
    end
end

[mean(Restimates) std(Restimates) length(Restimates) std(Restimates)/sqrt(length(Restimates))]
median(Restimates)
%mean(Restimates(Restimates<0.09))

%% Try reprocessing with a different erosion / dilation radius
erodeum=0.2; dilum =0.2;

Restimates = [];
for cc = 1:length(allcells)        
    CL = allcells(cc);   
    CL.reprocessROIs(erodeum,dilum);
    for sc = 1:length(CL.ROIgroups)       
        CL.ROIgroups(sc).Restimate = getRestimates(CL,CL.ROIgroups(sc),struct('mintubelen',5));
        Restimates(end+1) = CL.ROIgroups(sc).Restimate;
    end
end
Restimates = Restimates(~isnan(Restimates));
Restimates(Restimates>0.1) = [];

[mean(Restimates) std(Restimates) length(Restimates) std(Restimates)/sqrt(length(Restimates))]
median(Restimates)
%% Compare multiple datasets
load('200807_COS7_RTN4_KO_2G3_SNAP_KDEL_505_Sec61_Halo_TMR.mat')
erodeum=0.2; dilum =0.2;

RestimatesRTN4 = [];
Restimates1 = [];
for cc = 1:length(allcells)    
    CL = allcells(cc);
    CL.reprocessROIs(erodeum,dilum);
    for sc = 1:length(CL.ROIgroups)       
        CL.ROIgroups(sc).Restimate = getRestimates(CL,CL.ROIgroups(sc),struct('mintubelen',5));
        Restimates1(end+1) = CL.ROIgroups(sc).Restimate;
    end    
end
RestimatesRTN4 = [RestimatesRTN4 Restimates1];

load('../results/210106_COS7_RTN4_KO_Sec61_Halo_OG_ER_mCherry.mat')
Restimates2 = [];
for cc = 1:length(allcells)    
    CL = allcells(cc);    
    CL.reprocessROIs(erodeum,dilum);
    for sc = 1:length(CL.ROIgroups)       
        CL.ROIgroups(sc).Restimate = getRestimates(CL,CL.ROIgroups(sc),struct('mintubelen',5));
        Restimates2(end+1) = CL.ROIgroups(sc).Restimate;
    end        
end
RestimatesRTN4 = [RestimatesRTN4 Restimates2];


load('../results/210106_COS7_RTN4_KO_Sec61_Halo_OG_ER_mCherry_Lightning.mat');
Restimates3 = [];
for cc = 1:length(allcells)    
    CL = allcells(cc);    
    CL.reprocessROIs(erodeum,dilum);
    for sc = 1:length(CL.ROIgroups)       
        CL.ROIgroups(sc).Restimate = getRestimates(CL,CL.ROIgroups(sc),struct('mintubelen',5));
        Restimates3(end+1) = CL.ROIgroups(sc).Restimate;
    end        
end
RestimatesRTN4 = [RestimatesRTN4 Restimates3];

%%
RestimatesWT = [];
erodeum=0.2; dilum =0.2;

load('../results/210106_COS7_WT_Sec61_Halo_OG_ERmCherry.mat')
Restimates1 = [];
for cc = 1:length(allcells)    
    CL = allcells(cc);
    CL.reprocessROIs(erodeum,dilum);
    for sc = 1:length(CL.ROIgroups)       
        CL.ROIgroups(sc).Restimate = getRestimates(CL,CL.ROIgroups(sc),struct('mintubelen',5));
        Restimates1(end+1) = CL.ROIgroups(sc).Restimate;
    end            
end
RestimatesWT = [RestimatesWT Restimates1];

%
load('../results/200807_COS7_WT_SNAP_KDEL_505_Sec61_Halo_TMR.mat')
Restimates2 = [];
for cc = 1:length(allcells)    
    CL = allcells(cc);
    
    % exclude low resolution cells
    lowres = contains(CL.Name,'_1');
    if (lowres); continue; end
    CL.reprocessROIs(erodeum,dilum);
    for sc = 1:length(CL.ROIgroups)       
        CL.ROIgroups(sc).Restimate = getRestimates(CL,CL.ROIgroups(sc),struct('mintubelen',5));
        Restimates2(end+1) = CL.ROIgroups(sc).Restimate;
    end            
end
RestimatesWT = [RestimatesWT Restimates2];

load('../results/210106_COS7_WT_Sec61_Halo_OG_ER_mCherry_Lightning.mat')
Restimates3 = [];
for cc = 1:length(allcells)    
    CL = allcells(cc);
    
    % exclude low resolution cells   
    CL.reprocessROIs(erodeum,dilum);
    for sc = 1:length(CL.ROIgroups)       
        CL.ROIgroups(sc).Restimate = getRestimates(CL,CL.ROIgroups(sc),struct('mintubelen',5));
        Restimates3(end+1) = CL.ROIgroups(sc).Restimate;
    end            
end
RestimatesWT = [RestimatesWT Restimates3];

%%
[nanmedian(RestimatesWT) nanmedian(RestimatesRTN4)]


% --------------------------
%% get standard error from bootstrapping
ntrial = 500;
for bc = 1:ntrial
    valsWT = datasample(RestimatesWT,length(RestimatesWT),'Replace',true);
    valsRTN4 = datasample(RestimatesRTN4,length(RestimatesRTN4),'Replace',true);
    
    bootstrapWT(bc) = nanmedian(valsWT);
    bootstrapRTN4(bc) = nanmedian(valsRTN4);
end

disp('medians')
[median(RestimatesWT) median(RestimatesRTN4)]
disp('means of bootstrap')
[mean(bootstrapWT) mean(bootstrapRTN4)]
disp('std of bootstrap')
[std(bootstrapWT) std(bootstrapRTN4)]