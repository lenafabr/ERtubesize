dirname = '/data/proj/ERtransport/Edward081020/200807_COS7_RTN4 KO_lumern and membrane/extractedfiles/';
%fglob = '200807_COS7_WT_SNAP_KDEL_505_Sec61_Halo_TMR_Series*.tif'
fglob = '200807_COS7_RTN4_KO_2G3_SNAP_KDEL_505_Sec61_Halo_TMR*.tif';
files = dir([dirname fglob])

% offset for dilation and erosion, to deal with resolution issues
resoffset = 0.3;

%%
cellct = 0;
%%
%cellct = 0;
for fc = 5:length(files)
    fname = files(fc).name;
    [filepath,name,ext] = fileparts(fname);
    
    CL = CellObjTubeSheet(name);
    memframe = 2; lumframe = 1;
    CL.loadCellData(dirname,fname,memframe,lumframe);        
    
    % show membrane and luminal images
    figure(2)
    subplot(1,2,1)
    imshow(CL.imgmem,[])
    title('membrane')
    subplot(1,2,2)
    imshow(CL.imglum,[])
    title('luminal')
    
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
        sheetROI = CL.ROIgroups(sheetcount).sheetROI;
        tubeROIs = CL.ROIgroups(sheetcount).tubeROIs;
        Fsheet = sum(sum(sheetROI.erodemask.*double(CL.imgmem-CL.bgROI.avg)));
        Asheet = nnz(sheetROI.erodemask(:));
        
        Ftubes = zeros(length(tubeROIs),1);
        Ltubes = Ftubes;
        for tc = 1:length(tubeROIs)
            tube = tubeROIs(tc);
            Ftubes(tc) = sum(sum(tube.dilmask.*double(CL.imgmem-CL.bgROI.avg)));
            Ltubes(tc) = tube.L+2*tube.dilum*CL.pxperum;
        end
        
        totFtube = sum(Ftubes);
        totLtube = sum(Ltubes);
        
        Rtube = totFtube/Fsheet*Asheet/totLtube/pi/CL.pxperum;
        
        CL.ROIgroups(sheetcount).Restimate = Rtube;
        disp(sprintf('Group %d: Restimate %f', sheetcount,Rtube))
        
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
%save('../results/200807_COS7_WT_SNAP_KDEL_505_Sec61_Halo_TMR.mat')
save('200807_COS7_RTN4_KO_2G3_SNAP_KDEL_505_Sec61_Halo_TMR.mat')


%% get averages
Restimates = [];
for cc = 1:length(allcells)    
    CL = allcells(cc);
    
    % exclude low resolution cells
    lowres = contains(CL.Name,'_1');
    if (lowres); continue; end
    
    for sc = 1:length(CL.ROIgroups)        
        Restimates(end+1) = CL.ROIgroups(sc).Restimate;
    end
end

[mean(Restimates) std(Restimates) length(Restimates) median(Restimates)]
%mean(Restimates(Restimates<0.09))