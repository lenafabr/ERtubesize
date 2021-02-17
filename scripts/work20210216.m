% play around with ROI region
%load('../results/210106_COS7_WT_Sec61_Halo_OG_ERmCherry.mat')
%load('../results/200807_COS7_WT_SNAP_KDEL_505_Sec61_Halo_TMR.mat')
%load('../results/400ngCtrl_FITC_Calnexin_membrane.mat')
load('../results/210106_COS7_RTN4_KO_Sec61_Halo_OG_ER_mCherry.mat')
dilum=0.6; erodeum=0.3; endum=0.3
Restimates = [];
for cc = 1:length(allcells)    
    cc
    CL = allcells(cc);
    CL.reprocessROIsInterp(erodeum,dilum,endum);
    for sc = 1:length(CL.ROIgroups)       
        CL.ROIgroups(sc).Restimate = getRestimatesInterp(CL,CL.ROIgroups(sc),struct('mintubelen',5));
        if (~isnan(CL.ROIgroups(sc).Restimate))
            Restimates(end+1) = CL.ROIgroups(sc).Restimate;
        end
    end            
end

%% look at signals from particular tubule
%%
tube = CL.ROIgroups(1).tubeROIs(1);
%%
gaussfun = @(c,x) c(1)*exp(-(x-c(2)).^2./(2*c(3)^2)) + c(4);
%[tube] = processTubeInterp(CL,tube,dilum,endum)
linepos = linspace(-tube.dilum,tube.dilum,tube.nperp)*CL.pxperum;

for pc =2%:size(tube.interpvals,2)
    if (isnan(tube.intF(pc))); continue; end
    vals = tube.interpvals(:,pc-1);
    cfit = tube.gaussfit(:,pc-1);
    
    plot(linepos/CL.pxperum,vals,'.-',linepos/CL.pxperum,gaussfun(cfit,linepos),'k--')
    hold all
end
hold off

%%
CL.showROImasks(5)

%% look at signal perpendicular to tube
ROIgroup = CL.ROIgroups(1);
tube = ROIgroup.tubeROIs(2);

imshow(CL.imgmem,[])
hold all
path = tube.interppath;
plot(path(:,1),path(:,2),'.-')
hold off

%%
%CL = allcells(10);
ROIgroup = CL.ROIgroups(1);
tube = ROIgroup.tubeROIs(2);
path = tube.interppath;

ind = round(size(path,1)/2)

pt0 = path(ind,:);
pathv = path(ind,:)-path(ind-1,:);
pathv = pathv/norm(pathv);

pathperp = [pathv(2) -pathv(1)];
linepos = linspace(-8,8,10);
perppts = pt0 + linepos'*pathperp;

figure(1)
imshow(CL.imgmem,[])
hold all
path = tube.interppath;
plot(path(:,1),path(:,2),'.-')
plot(perppts(:,1),perppts(:,2),'.-')
hold off


[imgx imgy] = meshgrid(1:CL.ImgSize(1));
vals = interp2(imgx,imgy,double(CL.imgmem),perppts(:,1),perppts(:,2),'linear');

figure(2)
plot(linepos/CL.pxperum,vals,'.-')

%%
ROIgroup = CL.ROIgroups(1);
tube = ROIgroup.tubeROIs(3);
path = tube.interppath;

% estimate tube fluorescence from perpendicular projections
figure(1)
imshow(CL.imgmem,[])
hold all
path = tube.interppath;
plot(path(:,1),path(:,2),'r.-')

nperp = 20;
sumvals = zeros(nperp,1); perpct = 0;

linepos = linspace(-1,1,nperp)*CL.pxperum;

for pc = 15%2:size(path,1)-1
    pc
    pt0 = path(pc,:);
    pathv = path(pc,:)-path(pc-1,:);
    pathv = pathv/norm(pathv);

    pathperp = [pathv(2) -pathv(1)];    
    perppts = pt0 + linepos'*pathperp;
    
    [imgx imgy] = meshgrid(1:CL.ImgSize(1));
    vals = interp2(imgx,imgy,double(CL.imgmem),perppts(:,1),perppts(:,2),'linear');
    
    figure(1)
    plot(perppts(:,1),perppts(:,2),'.-')
    
    figure(2)
    plot(linepos/CL.pxperum,vals,'.-')
    hold all
    
    sumvals = sumvals + vals;
    perpct = perpct+1;
end

figure(1); hold off
figure(2); hold off

%%
figure(2)
hold all
plot(linepos/CL.pxperum,sumvals/perpct,'k.-','LineWidth',2)
hold off

%% integrate fluorescence and approximate radius
avgvals = sumvals/perpct;
dx = linepos(2)-linepos(1);
avgends = (avgvals(1)+avgvals(end))/2;
avgvals = avgvals - avgends;
intF = (sum(avgvals) - 0.5*(avgvals(1)+avgvals(end)))*dx;

sheetF = sum(sum(ROIgroup.sheetROI.erodemask.*double(CL.imgmem - avgends)))
sheetA = nnz(ROIgroup.sheetROI.erodemask);

intF/sheetF*sheetA/pi/CL.pxperum

%% fit gaussian to each tubule point

tube = CL.ROIgroups(1).tubeROIs(2);
dilum = 0.8; endum = 0.2;
[tube] = processTubeInterp(CL,tube,dilum,endum)

%
%tube = CL.ROIgroups(6).tubeROIs(1);
%
linepos = linspace(-tube.dilum,tube.dilum,20)*CL.pxperum;

gaussfun = @(c,x) c(1)*exp(-(x-c(2)).^2./(2*c(3)^2)) + c(4);
for pc = 1:size(tube.interpvals,2)
    plot(linepos/CL.pxperum,tube.interpvals(:,pc),...
        linepos/CL.pxperum,gaussfun(tube.gaussfit(:,pc),linepos),'k--')
    hold all
end
hold off

%%
tube = CL.ROIgroups(5).tubeROIs(1);

gaussfun = @(c,x) c(1)*exp(-(x-c(2)).^2./(2*c(3)^2)) + c(4);
%[tube] = processTubeInterp(CL,tube,dilum,endum)
linepos = linspace(-tube.dilum,tube.dilum,tube.nperp)*CL.pxperum;

for pc = 11%:size(tube.interpvals,2)
    vals = tube.interpvals(:,pc);
    cfit = tube.gaussfit(:,pc);
    
    plot(linepos/CL.pxperum,vals,'.-',linepos/CL.pxperum,gaussfun(cfit,linepos),'k--')
    hold all
end
hold off