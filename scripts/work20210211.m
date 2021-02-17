%dirname = '/data/proj/ERtransport/210105_COS7_RTN4 KO_tube narrowing_STELLARIS8/extractedimages/';
%filename = 'Series009.tif';
dirname = '/data/proj/ERtransport/Edward081020/200807_COS7_RTN4 KO_lumern and membrane/extractedfiles/';
filename = '200807_COS7_WT_SNAP_KDEL_505_Sec61_Halo_TMR_Series003_2.tif'


img = imread([dirname filename],2);
imglum = imread([dirname filename],1);
info = imfinfo([dirname filename]);
pxperum = info.XResolution;
%%
subplot(1,2,1)
imshow(img,[])
title('membrane')
subplot(1,2,2)
imshow(imglum,[])
title('luminal')

% offset in um to deal with resolution
resoffset = 0.3;

%% compute background
subplot(1,2,1)
bgroi = drawpolygon()
%
bgmask = createMask(bgroi);
% get bg signal per px
bg = sum(sum(bgmask.*double(img)))/nnz(bgmask(:))

%% manually select sheet region
subplot(1,2,1)
sheetroi = drawpolygon()
%%
sheetbound = sheetroi.Position;
%%
% get sheet mask and erode
sheetmask = createMask(sheetroi);

subplot(1,2,1)
imshowpair(img,sheetmask)

%% erode sheet by some amount on either side
px = resoffset*pxperum;
se = strel('disk',floor(px))

sheeterode = imerode(sheetmask,se);
subplot(1,2,1)
imshowpair(img,sheeterode)

%%
% draw tubules
subplot(1,2,1)
imshow(img,[])
roitube = drawpolyline()
%%

tubepath = roitube.Position;


%% create mask along tube path

% interpolate path to denser points
[param,lens] = arclenparam(tubepath');
npt = ceil(param(end));
newparam = linspace(2*resoffset*pxperum,max(param)-2*resoffset*pxperum,npt);
%newparam = linspace(0,max(param));
interppath = interp1(param,tubepath,newparam,'spline');
finparam = arclenparam(interppath');
    
% plot(tubepath(:,1),tubepath(:,2),'.-')
% hold all
% plot(interppath(:,1),interppath(:,2),'.-')
% hold off
%
% imshow(img,[])
% hold all
% plot(interppath(:,1),interppath(:,2),'r.-')
% hold off
%
tubepx = round(interppath);
ind= sub2ind(size(img),tubepx(:,2),tubepx(:,1))
tubemask = zeros(size(img));
tubemask(ind)=1;

px = resoffset*pxperum;
se = strel('disk',floor(px))
tubedilate = imdilate(tubemask,se);

subplot(1,2,1)
imshowpair(img,tubedilate)
hold all
plot(interppath(:,1),interppath(:,2),'c.-')
hold off

%% get brightness in sheets and tubules

Fsheet = sum(sum(sheeterode.*double(img-bg)));
Asheet = nnz(sheeterode(:));
%Fsheet = sum(sum(sheetmask.*double(img)));
%Asheet = nnz(sheetmask(:));

Ftube = sum(sum(tubedilate.*double(img-bg)));
Ltube = finparam(end);


Rtube = Ftube/Fsheet*Asheet/(Ltube+2*resoffset*pxperum)/pi/pxperum


%% load in a MIP from stack from laura
dirname = '/data/proj/ERsorting/Cell 12_07.16.2019_12/';
fname = 'FITC_Cell12_Calnexin_membrane.tif';

info = imfinfo([dirname fname]);
nframe = length(info);
sz = info(1).Width;
pxperum = info(1).XResolution;
imgs = zeros(sz,sz,nframe);
for fc = 1:nframe
    imgs(:,:,fc) = imread([dirname fname],fc);
end

%%
mip = max(imgs,[],3);

imshow(mip,[])