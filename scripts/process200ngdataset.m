% ctrl cells
%dirname = '/data/proj/ERsorting/200 ng ER Shaping Proteins/01 COS7 Control/200 ng mCherry N1 200 ng mEmerald Calnexin 200 ng BFP KDEL/';
%fglob = 'Cell*2019*/*alnexin_membrane.tif';
% Rtn4OE cells
dirname = '/data/proj/ERsorting/200 ng ER Shaping Proteins/02 COS7 Rtn4a/200 ng mCherry Rtn4a 200 ng mEmerald Calnexin 200 ng BFP KDEL/';
fglob = 'Cell*2019*/*alnexin_membrane.tif';

files = dir([dirname fglob])

% offset for dilation and erosion, to deal with resolution issues
resoffset = 0.2;

%%
cellct = 0;
%%
%cellct = 0;
for fc = 18:length(files)
    fname = files(fc).name
    [filepath,name,ext] = fileparts(fname);
    
    CL = CellObjTubeSheet(name);
    
    %%
    CL.loadCellStackMaxI([files(fc).folder '/'], fname)
   
    %% adjust image
    img = CL.imgmem - min(CL.imgmem(:));
    img = img./max(img(:));
    img2 = imadjust(img,[0,0.3],[0,1]);
    CL.imgmem = img2;
    
    % show membrane and luminal images   
    figure(1)
    imshow(CL.imgmem,[])
    title('membrane')
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
        %imshowpair(CL.imgmem,sheetROI.erodemask)
        %imshowpair(img2,sheetROI.erodemask)
        
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
        CL.showROImasks()
        
                
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
%save('../results/200ngRTN4OE_FITC_Calnexin_membrane.mat', '-v7.3');
save('../results/200ngRTN4OE_Calnexin_membrane.mat', '-v7.3');
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

[mean(Restimates) std(Restimates) length(Restimates) std(Restimates)/sqrt(length(Restimates))]
median(Restimates)

%% save without images
% for cc = 1:length(allcells)
%     allcells(cc).imgmem = [];
% end
% save('../results/400ngCtrl_FITC_Calnexin_membrane_noimg.mat', 'allcells','dirname','fglob','files','dilum','erodeum');

%% load Laura's data

load('../results/400ngCtrl_b_FITC_Calnexin_membrane.mat')
%%
erodeum=0.2; dilum =0.2;
RestimatesWTLaura = [];
for cc = 1:length(allcells)    
    CL = allcells(cc);
    CL.reprocessROIs(erodeum,dilum);
    options = struct('mintubelen',5,'maxtubedist',4);    
    for sc = 1:length(CL.ROIgroups)       
        CL.ROIgroups(sc).Restimate = getRestimates(CL,CL.ROIgroups(sc),options);
        RestimatesWTLaura(end+1) = CL.ROIgroups(sc).Restimate;
    end    
end
%
RestimatesWTLaura = RestimatesWTLaura(~isnan(RestimatesWTLaura));

[mean(RestimatesWTLaura) median(RestimatesWTLaura) length(RestimatesWTLaura)]
%%
load('../results/400ngRTN4OE_b_FITC_Calnexin_membrane.mat')
erodeum=0.2; dilum =0.2;
RestimatesRTN4OE = [];
for cc = 1:length(allcells)    
    CL = allcells(cc);
    CL.reprocessROIs(erodeum,dilum);
    for sc = 1:length(CL.ROIgroups)       
        CL.ROIgroups(sc).Restimate = getRestimates(CL,CL.ROIgroups(sc),struct('mintubelen',5));
        RestimatesRTN4OE(end+1) = CL.ROIgroups(sc).Restimate;
    end    
end
%%
load('../results/200ngCtrl_Calnexin_membrane.mat')
%
erodeum=0.2; dilum =0.2;
Restimates200WTLaura = [];
options = struct('mintubelen',3,'maxtubedist',inf);    
for cc = 1:length(allcells)    
    CL = allcells(cc);
    CL.reprocessROIs(erodeum,dilum);
    for sc = 1:length(CL.ROIgroups)              
        CL.ROIgroups(sc).Restimate = getRestimates(CL,CL.ROIgroups(sc),options);
        Restimates200WTLaura(end+1) = CL.ROIgroups(sc).Restimate;
    end    
end
Restimates200WTLaura = Restimates200WTLaura(~isnan(Restimates200WTLaura));

[mean(Restimates200WTLaura) median(Restimates200WTLaura) length(Restimates200WTLaura)]
%%
load('../results/200ngRTN4OE_Calnexin_membrane.mat')
erodeum=0.2; dilum =0.2;
Restimates200RTN4OE = [];
options = struct('mintubelen',3,'maxtubedist',inf);   
for cc = 1:length(allcells)    
    CL = allcells(cc);
    CL.reprocessROIs(erodeum,dilum);
    for sc = 1:length(CL.ROIgroups)       
        CL.ROIgroups(sc).Restimate = getRestimates(CL,CL.ROIgroups(sc),options);
        Restimates200RTN4OE(end+1) = CL.ROIgroups(sc).Restimate;
    end    
end
Restimates200RTN4OE = Restimates200RTN4OE(~isnan(Restimates200RTN4OE));

[mean(Restimates200RTN4OE) median(Restimates200RTN4OE) length(Restimates200RTN4OE)]



%% Compare with Edward's datasets
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

%%
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

%%
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
RestimatesRTN4sr = Restimates3;
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

%
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
RestimatesWTsr = Restimates3;
RestimatesWT = [RestimatesWT Restimates3];

%%
[nanmedian(RestimatesWT) nanmedian(RestimatesRTN4OE)]

%% Compare all results with statistical tests

[hWT,pWT] = kstest2(RestimatesWT, RestimatesWTLaura)
[hRTN4OE, pRTN4OE] = kstest2(RestimatesRTN4OE, RestimatesWTLaura)
[hRTN4KO, pRTN4KO] = kstest2(RestimatesRTN4, RestimatesWT)
[hSR, pSR] = kstest2(RestimatesWTsr, RestimatesWT(1:end-length(RestimatesWTsr)))

%%
display({'Edward WT', 'Laura WT', 'Edward RTN4KO', 'Laura RTN4OE'})

[nanmean(RestimatesWT) nanmean(RestimatesWTLaura) nanmean(RestimatesRTN4) nanmean(RestimatesRTN4OE)]

ind = find(~isnan(RestimatesWT)); nWT = length(ind);
steWT = std(RestimatesWT(ind))/sqrt(nWT);
ind = find(~isnan(RestimatesWTLaura)); nWTLaura = length(ind);
steWTLaura = std(RestimatesWTLaura(ind))/sqrt(length(ind));
ind = find(~isnan(RestimatesRTN4)); nRTN4 = length(ind);
steRTN4 = std(RestimatesRTN4(ind))/sqrt(length(ind));
ind = find(~isnan(RestimatesRTN4OE)); nRTN4OE = length(ind);
steRTN4OE = std(RestimatesRTN4OE(ind))/sqrt(length(ind));

[steWT steWTLaura steRTN4 steRTN4OE]
[nWT nWTLaura nRTN4 nRTN4OE]


disp('super resolution WT and RTN4KO')
[nanmean(RestimatesWTsr) nanmean(RestimatesRTN4sr)]
ind = find(~isnan(RestimatesWTsr)); nWTsr = length(ind);
steWTsr = std(RestimatesWTsr(ind))/sqrt(length(ind));
ind = find(~isnan(RestimatesRTN4sr)); nRTN4sr = length(ind);
steRTN4sr = std(RestimatesRTN4sr(ind))/sqrt(length(ind));
[steWTsr steRTN4sr]
[nWTsr nRTN4sr]
% --------------------------
%% get standard error from bootstrapping
ntrial = 500;
for bc = 1:ntrial
    valsWT = datasample(RestimatesWT,length(RestimatesWT),'Replace',true);
    valsRTN4 = datasample(RestimatesRTN4,length(RestimatesRTN4),'Replace',true);
    
    bootstrapWT(bc) = nanmean(valsWT);
    bootstrapRTN4(bc) = nanmean(valsRTN4);
end

disp('medians')
[median(RestimatesWT) nanmedian(RestimatesRTN4)]
disp('means of bootstrap')
[mean(bootstrapWT) mean(bootstrapRTN4)]
disp('std of bootstrap')
[std(bootstrapWT) std(bootstrapRTN4)]

% ----------
%% Try keeping only sufficiently nearby tubules
