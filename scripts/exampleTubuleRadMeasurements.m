%dirname = ['~/UCSD/data/Avezov/Tasuko/220421_COS7_Halo-Sec61b-TMR_for tubule width_live and fixed/'];
dirname = ['/data/proj/ERtransport/Tasuko20220421_tubewidth/220421_COS7_Halo-Sec61b-TMR_for tubule width_live and fixed/'];
%
fglob = ['220421_COS7_Halo_Sec61b_TMR_fortubulewidth_live*.tif'];


files = dir([dirname fglob])

%% offset for dilation and erosion, to deal with resolution issues
resoffset = 0.2;
cellct = 0;

%%
for fc = 3:length(files)
    fname = files(fc).name
    [filepath,name,ext] = fileparts(fname);
    
    CL = CellObjTubeSheet(name);
    memframe = 1; lumframe = 1;
    CL.loadCellData(dirname,fname,memframe,lumframe);        
    CL.imgmem = double(CL.imgmem)/max(double(CL.imgmem(:)));
    % brighten the image by maxing out top percentage of pixels
    % to not adjust at all, put [0,1],[0,1] for the arguments
    imgmemadj = imadjust(CL.imgmem,[0,0.75],[0,1]);
    
    % show membrane and luminal images
    figure(1)    
    imshow(imgmemadj)
    title(sprintf('Cell %d: %s',fc,fname))
    
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
        imshowpair(imgmemadj,sheetROI.erodemask)
            
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
%save(['../results/220421_COS7_Halo_Sec61b_TMR_fortubulewidth_live.mat'])
save(['../results/example.mat'])

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
