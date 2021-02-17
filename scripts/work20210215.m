% play around with ROI region
load('../results/210106_COS7_WT_Sec61_Halo_OG_ERmCherry.mat')
dilum=0.3; erodeum=0.3;
Restimates = [];
for cc = 1:length(allcells)    
    CL = allcells(cc);
    CL.reprocessROIs(erodeum,dilum);
    for sc = 1:length(CL.ROIgroups)       
        CL.ROIgroups(sc).Restimate = getRestimates(CL,CL.ROIgroups(sc),struct('mintubelen',5));
        Restimates(end+1) = CL.ROIgroups(sc).Restimate;
    end            
end
%%
CL.showROImasks(5)

%% look at signal perpendicular to tube
ROIgroup = CL.ROIgroups(5);
tube = ROIgroup.tubeROIs(3);

imshow(CL.imgmem,[])
hold all
path = tube.interppath;
plot(path(:,1),path(:,2),'.-')
hold off

%%
%CL = allcells(10);
ROIgroup = CL.ROIgroups(5);
tube = ROIgroup.tubeROIs(4);
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
[param] = arclenparam(tubeROI.path');
            npt = max(2,ceil(param(end)-2*dilum*CL.pxperum));
            %newparam = linspace(2*dilum*CL.pxperum,max(param)-2*dilum*CL.pxperum,npt);
            newparam = linspace(dilum*CL.pxperum,max(param)-dilum*CL.pxperum,npt);
            %newparam = linspace(0,max(param),npt);
            tubeROI.interppath = interp1(param,tubeROI.path,newparam,'spline');

%%
ROIgroup = CL.ROIgroups(2);

sheetROI = ROIgroup.sheetROI;

tubeROIs = ROIgroup.tubeROIs;

%%
tubeROI = tubeROIs(2);
imshowpair(CL.imgmem,tube.mask)
hold all
path = tubeROI.interppath;
plot(path(:,1),path(:,2),'y.-')
hold off

%%
dilum = 0.3;
% interpolate path to denser points
[param] = arclenparam(tubeROI.path');
npt = ceil(param(end));
newparam = linspace(0,max(param),npt);
tubeROI.interppath = interp1(param,tubeROI.path,newparam,'spline');
finparam = arclenparam(tubeROI.interppath');
tubeROI.L = finparam(end);
tubeROI.dilum = dilum;


% dilate and make mask
tubepx = round(tubeROI.interppath);
ind= sub2ind(CL.ImgSize,tubepx(:,2),tubepx(:,1));
tubeROI.mask = zeros(CL.ImgSize);
tubeROI.mask(ind)=1;

px = dilum*CL.pxperum;
se = strel('disk',floor(px));
tubeROI.dilmask = imdilate(tubeROI.mask,se);

imshowpair(CL.imgmem,tubeROI.dilmask)
hold all
path = tubeROI.interppath;
plot(path(:,1),path(:,2),'y.-')
hold off

%% project mask to curve and get rid of any beyond ends
maskind = find(tubeROI.dilmask);
[maskx, masky] = ind2sub(CL.ImgSize,maskind);
[xy,distance,t_a] = distance2curve(tubeROI.path,[masky maskx]);
dropind = find(t_a<=0 | t_a>=1);
newmask = tubeROI.dilmask;
newmask(maskind(dropind))=0;

imshowpair(CL.imgmem,newmask)
hold all
path = tubeROI.interppath;
plot(path(:,1),path(:,2),'y.-')
hold off

%% look at signal perpendicular to tube
