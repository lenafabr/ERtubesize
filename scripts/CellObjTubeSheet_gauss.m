% class definition for a cell object, to be used
% in calculating tube width

classdef CellObjTubeSheet < handle
    
    properties
        % Cell Name
        Name = 'default';
        
        % image size
        ImgSize = [];
        
        % directory where image data is stored
        DirName = '';
        
        % filename:
        filename = ''; 
                
        % spatial resolution (px/um)
        pxperum = 0;
        
        % single snapshot image of ER membrane and lumen in this cell
        % uses up memory but not too much since it's just one snapshot
        imgmem = NaN;
        imglum = NaN;
        
        % array of info about sheet and corresponding tube ROIs;
        ROIgroups = [];         
        
        % background ROI info, including bg (background signal) and mask        
        bgROI = NaN;                  
        
    end
    
    methods
        function CL = CellObjTubeSheet(name)
            % create a cell object with the given name
            CL.Name = name;
            
            % set up empty ROIs structure
            % each ROIgroup includes a single sheet ROI and multiple tube ROIs 
            CL.ROIgroups = struct('sheetROI',{},'tubeROIs',{},'Restimate',{})
            %= struct('ID',{},'bound',{},'mask',{},'erodemask',{},'erodeum',{});
            %CL.tubeROIs = struct('ID',{},'mask',{},'dilmask',{},'L',{},'dilum',{},'sheetid',{},'path',{},'interppath',{});            
        end
        
%         function clearTubeROIs(CL)
%             CL.tubeROIs = struct('ID',{},'mask',{},'dilmask',{},'L',{},'dilum',{},'sheetid',{},'path',{},'interppath',{});            
%         end
        function loadCellData(CL,dirname,filename,memframe,lumframe)
            % load in image data
            CL.DirName = dirname;
            CL.filename = filename;
            
            CL.imgmem = imread([dirname filename],memframe);
            CL.imglum = imread([dirname filename],lumframe);
            
            info = imfinfo([dirname filename]);
            CL.pxperum = info.XResolution;
            
            CL.ImgSize = size(CL.imgmem);            
        end
        
        function loadCellStackMaxI(CL,dirname, filename)
            % wrote a stack of images and get max intensity projection
            CL.DirName = dirname;
            CL.filename = filename;
            
            info = imfinfo([dirname filename]);
            CL.pxperum = info(1).XResolution;
            imgs = zeros(info(1).Width, info(1).Width,length(info));
            for ic = 1:length(info)
                imgs(:,:,ic) = imread([dirname filename],ic);
            end
            
            CL.imgmem = max(imgs,[],3);
            
            CL.ImgSize = size(CL.imgmem);
        end
        
        function getBackground(CL,fig)
            % pick out a background ROI and calculate bg signal per px
            
            if (exist('fig','var'))
                figure(fig)
            else
                figure
            end
            
            disp('pick background region')
            bgroi = drawpolygon()
            input('hit enter when done')
            
            CL.bgROI = struct();
            CL.bgROI.bound = bgroi.Position;            
            CL.bgROI.mask = createMask(bgroi);
            
            % get bg signal per px
            CL.bgROI.avg = sum(sum(CL.bgROI.mask.*double(CL.imgmem)))/nnz(CL.bgROI.mask(:))            
        end
        
        function sheetROI = getSheet(CL,erodeum,fig)
            % pick out a sheet ROI
            %  provide um for erosion 
            
            if (exist('fig','var'))
                figure(fig)
            else
                figure
                imshow(CL.imgmem,[])
            end            
            
            disp('pick sheet region')
            sheetroi = drawpolygon()            
            input('Hit enter when done.\n')
                      
            sheetROI = struct();
            sheetROI.bound = sheetroi.Position;
            sheetROI.mask = createMask(sheetroi);
            sheetROI.erodeum = erodeum;

            px = erodeum*CL.pxperum;
            se = strel('disk',floor(px));

            sheetROI.erodemask = imerode(sheetROI.mask,se);
            sheetROI.tubeids = [];        
        end
        
        
        function tubeROI = getTube(CL,dilum,fig)
            % pick out a tube ROI, attached to a particular sheet
            % dilum give dilation um
            
            if (exist('fig','var'))
                figure(fig)
            else
                figure
                imshow(CL.imgmem,[])
            end            
            
            disp('pick tube region')
            roitube = drawpolyline()
            input('Hit enter when done.\n')
            tubeROI = struct();
            
            tubeROI.path = roitube.Position;
            tubeROI = CL.processTube(tubeROI,dilum);                                          
        end
        
        function tubeROI= processTube(CL,tube,dilum)
            % process a tube ROI object, calculating dilation amsk
            % path must already be set
            
            tubeROI = tube;
            
            % interpolate path to denser points
            [param] = arclenparam(tubeROI.path');
            npt = max(2,ceil(param(end)-2*dilum*CL.pxperum));
            %newparam = linspace(2*dilum*CL.pxperum,max(param)-2*dilum*CL.pxperum,npt);
            newparam = linspace(dilum*CL.pxperum,max(param)-dilum*CL.pxperum,npt);
            %newparam = linspace(0,max(param),npt);
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
            
            maskind = find(tubeROI.dilmask);
            [maskx, masky] = ind2sub(CL.ImgSize,maskind);
            [xy,distance,t_a] = distance2curve(tubeROI.interppath,[masky maskx]);
            dropind = find(t_a<=0 | t_a>=1);
            newmask = tubeROI.dilmask;
            newmask(maskind(dropind))=0;
            tubeROI.dilmask = newmask;
            
        end
        
        function [tubeROI] = processTubeInterp(CL,tube,dilum,endum)
            % process a tube ROI object, calculating dilation amsk
            % path must already be set
            % use subpixel interpolation perpendicular to tubule path
            % endum = how much to cut off end
            
            tubeROI = tube;
            
            % interpolate path to denser points
            [param] = arclenparam(tubeROI.path');
            npt = max(2,ceil(param(end)-2*endum*CL.pxperum));
            %newparam = linspace(2*dilum*CL.pxperum,max(param)-2*dilum*CL.pxperum,npt);
            newparam = linspace(endum*CL.pxperum,max(param)-endum*CL.pxperum,npt);
            %newparam = linspace(0,max(param),npt);
            tubeROI.interppath = interp1(param,tubeROI.path,newparam,'spline');
            finparam = arclenparam(tubeROI.interppath');
            tubeROI.L = finparam(end);
            tubeROI.dilum = dilum;
            
            
            % make mask
            tubepx = round(tubeROI.interppath);
            ind= sub2ind(CL.ImgSize,tubepx(:,2),tubepx(:,1));
            tubeROI.mask = zeros(CL.ImgSize);
            tubeROI.mask(ind)=1;
            
            
            % interpolated dilation
            dilmask = zeros(size(tubeROI.mask));
            nperp = 20;
            sumvals = zeros(nperp,1); perpct = 0;
            
            linepos = linspace(-dilum,dilum,nperp)*CL.pxperum;            
            path = tubeROI.interppath;
            [imgx imgy] = meshgrid(1:CL.ImgSize(1));
            
            tubeROI.interpvals = zeros(nperp,size(path,1)-1);
            tubeROI.gaussfit = zeros(4,size(path,1)-1);
            
            startc = zeros(1,4);
            intF = zeros(size(path,1)-1,1);
            for pc = 2:size(path,1)
                pt0 = path(pc,:);
                pathv = path(pc,:)-path(pc-1,:);
                normstep(pc-1) = norm(pathv);
                pathv = pathv/normstep(pc-1);
                
                pathperp = [pathv(2) -pathv(1)];
                perppts = pt0 + linepos'*pathperp;
               
                perppts = max(perppts,1);
                perppts  = min(perppts,CL.ImgSize(1));
                
                vals = interp2(imgx,imgy,double(CL.imgmem),perppts(:,1),perppts(:,2),'linear');
                
                tubeROI.interpvals(:,pc-1) = vals;
                
                % find minimum on either side of zero and do not let it
                % grow again
                [~,m1] = min(vals(1:nperp/2));
                [~,m2] = min(vals(nperp/2+1:end)); m2 = m2+nperp/2;               
                vals(1:m1-1) = vals(m1);
                vals(m2+1:end) = vals(m2);
                
                % fit to a gaussian
                gaussfun = @(c,x) c(1)*exp(-(x-c(2)).^2./(2*c(3)^2)) + c(4);        
                startc(1) = max(vals)-min(vals);
                startc(2) = 0;
                startc(3) = dilum/2;
                startc(4) = min(vals);
                [cfit,Res] = nlinfit(linepos,vals',gaussfun,startc);
                cfit(3) = abs(cfit(3));
                
                meanRes = mean(Res)./mean(vals); % fractional mean residual                                
                               
                tubeROI.gaussfit(:,pc-1) = cfit';
                tubeROI.nperp = nperp;                                
                localbg(pc) = cfit(4);
                
                % check if local min on either side of peak is too big
                % sign of bad fit or overlap with another feature
                % minimum should be less than half of the way from bottom
                % to top
                [~,ind] = min(abs(linepos-cfit(2)));
                diff1 = abs(min(vals(1:floor(ind))) - cfit(4));
                diff2 = abs(min(vals(ceil(ind):end)) - cfit(4));
                badfit = max(diff1,diff2)>cfit(1)*0.5;                                              
                
                if (badfit || abs(meanRes)>1 || abs(cfit(2))>linepos(end)/2)
                    %tubeROI.gaussfit(:,pc-1) = NaN*cfit';     
                    localbg(pc) = NaN;
                    intF(pc-1) = NaN;                    
                    continue
                end
                
                sumvals = sumvals + vals;
                perpct = perpct+1;
                
                % integrate fluorescence over +/- 2 sigma
                xperp = linspace(cfit(2)-cfit(3)*2,cfit(2)+cfit(3)*2,nperp);
                perppts = pt0 + xperp'*pathperp;
                perppts = max(perppts,1);
                perppts  = min(perppts,CL.ImgSize(1));
                % adjust point to center of gaussian fit
                % tubeROI.interppath(pc,:) = pt0 + cfit(2)*pathperp;
                valint = interp2(imgx,imgy,double(CL.imgmem - localbg(pc)),perppts(:,1),perppts(:,2),'linear');
                                
                ind = sub2ind(CL.ImgSize,floor(perppts(:,2)),floor(perppts(:,1)));
                dilmask(ind) = 1;                
                
                %valint = interp1(linepos,vals,xperp) - localbg(pc);
                intF(pc-1) = (sum(valint)-0.5*(valint(1)+valint(end)))*normstep(pc-1)*(xperp(2)-xperp(1));        
            end
            
            % estimate local background
            tubeROI.localbg = mean(localbg);
            
            tubeROI.dilmask = dilmask;
            tubeROI.intF = intF;
            tubeROI.localbg = mean(tubeROI.gaussfit(4,:));
            
            badind = find(isnan(intF));            
           % tubeROI.intF(badind) = [];
           % tubeROI.L = tube.L-sum(normstep(badind));
            tubeROI.normstep = normstep;
%            tubeROI.gaussfit(:,badind) = [];
%            avgvals = sumvals/perpct;
%             
%             dx = linepos(2)-linepos(1);
%             % estimate local background
%             avgends = (avgvals(1)+avgvals(end))/2;
%             avgvals = avgvals - avgends;
%             tubeROI.intF = (sum(avgvals) - 0.5*(avgvals(1)+avgvals(end)))*dx;
%             tubeROI.localbg = avgends;
        end
        
        function reprocessROIsInterp(CL,erodeum,dilum,endum)
            % use interpolation to reprocess all ROIs
            for rc = 1:length(CL.ROIgroups)
                sheetROI = CL.ROIgroups(rc).sheetROI;
                sheetROI.erodeum = erodeum;
                px = erodeum*CL.pxperum;
                se = strel('disk',floor(px));
                sheetROI.erodemask = imerode(sheetROI.mask,se);
               
                CL.ROIgroups(rc).sheetROI = sheetROI;
                tubeROIs = CL.ROIgroups(rc).tubeROIs;
                                
                
                for tc = 1:length(tubeROIs)
                    tubeROI = tubeROIs(tc);                    
                    tubeROI = CL.processTubeInterp(tubeROI,dilum,endum);
                    if (size(tubeROI.interppath,1) > 1)
                        newtubeROIs(tc) = tubeROI;
                    end                    
                end
                CL.ROIgroups(rc).tubeROIs = newtubeROIs;
                
                sheetROI.localbg = mean([CL.ROIgroups(rc).tubeROIs.localbg]);
                CL.ROIgroups(rc).sheetROI = sheetROI;
            end
        end
        
        function reprocessROIs(CL,erodeum,dilum)
            %% reprocess ROIs using a new value for the erosion and dilution
            % radius
            
            for rc = 1:length(CL.ROIgroups)
                sheetROI = CL.ROIgroups(rc).sheetROI;
                sheetROI.erodeum = erodeum;
                px = erodeum*CL.pxperum;
                se = strel('disk',floor(px));
                sheetROI.erodemask = imerode(sheetROI.mask,se);
               
                CL.ROIgroups(rc).sheetROI = sheetROI;
                
                for tc = 1:length(CL.ROIgroups(rc).tubeROIs)
                    tubeROI = CL.ROIgroups(rc).tubeROIs(tc);                    
                    tubeROI = CL.processTube(tubeROI,dilum);                                        
                    CL.ROIgroups(rc).tubeROIs(tc) = tubeROI;
                end
            end
        end
        
        function Rtube = getRestimates(CL,ROIgroup,options)
            % get radius estimate for a particular ROI group
            
            opt = struct();
            % minimal tube length in pixels
            opt.mintubelen = 0;
            
            if (exist('options','var'))
                opt = copyStruct(options,opt);
            end
            
            %% calculate radius estimate for this sheet and tubule group
            sheetROI = ROIgroup.sheetROI;
            tubeROIs = ROIgroup.tubeROIs;
            Fsheet = sum(sum(sheetROI.erodemask.*double(CL.imgmem-CL.bgROI.avg)));
            Asheet = nnz(sheetROI.erodemask(:));
            
            Ftubes = zeros(length(tubeROIs),1);
            Ltubes = Ftubes;
            for tc = 1:length(tubeROIs)
                tube = tubeROIs(tc);
                Ftubes(tc) = sum(sum(tube.dilmask.*double(CL.imgmem-CL.bgROI.avg)));
                Ltubes(tc) = tube.L;%+2*tube.dilum*CL.pxperum;
            end
            
            if (opt.mintubelen>0)
                goodind = find(Ltubes>opt.mintubelen);
            else
                goodind = 1:length(Ltubes);
            end
            
            totFtube = sum(Ftubes(goodind));
            totLtube = sum(Ltubes(goodind));
            
            Rtube = totFtube/Fsheet*Asheet/totLtube/pi/CL.pxperum;            
        end
        
        function Rtube = getRestimatesInterp(CL,ROIgroup,options)
            % get radius estimate for a particular ROI group
            % using interpolated tube fluorescence
            
            opt = struct();
            % minimal tube length in pixels
            opt.mintubelen = 0;
            
            if (exist('options','var'))
                opt = copyStruct(options,opt);
            end
            
            %% calculate radius estimate for this sheet and tubule group
            sheetROI = ROIgroup.sheetROI;
            tubeROIs = ROIgroup.tubeROIs;
            Fsheet = sum(sum(sheetROI.erodemask.*double(CL.imgmem-CL.bgROI.avg)));
            %Fsheet = sum(sum(sheetROI.erodemask.*double(CL.imgmem-sheetROI.localbg)));
            Asheet = nnz(sheetROI.erodemask(:));
            
            %Ltubes = [tubeROIs.L];
            for tc = 1:length(tubeROIs)
                goodind = ~isnan(tubeROIs(tc).intF);
                Ftubes(tc) = sum(tubeROIs(tc).intF(goodind));
                Ltubes(tc) = sum(tubeROIs(tc).normstep(goodind));
            end            
            
            if (opt.mintubelen>0)
                goodind = find(Ltubes>opt.mintubelen);
            else
                goodind = 1:length(Ltubes);
            end
            
            %totFL = sum(Ftubes(goodind).*Ltubes(goodind));
            totFL = sum(Ftubes(goodind));
            totLtube = sum(Ltubes(goodind));
            
            Rtube = totFL/Fsheet*Asheet/totLtube/pi/CL.pxperum;            
        end
        
        function showROImasks(CL,whichroi,showlum)
            
            if (whichroi<=0)
                whichroi = 1:length(CL.ROIgroups);
            end
            
            if (~exist('showlum','var'))
                showlum = false;
            end
            
            %% show all sheet and tube masks and paths
            totmask = zeros(CL.ImgSize);
            for sc= whichroi
                sheetROI = CL.ROIgroups(sc).sheetROI;
                tubeROIs = CL.ROIgroups(sc).tubeROIs;
                totmask = totmask + sheetROI.erodemask;
                paths = {};
                for tc = 1:length(tubeROIs)
                    totmask = totmask + tubeROIs(tc).dilmask;
                end
            end
            if (showlum)
                imshowpair(CL.imglum,totmask)
            else
                imshowpair(CL.imgmem,totmask)
            end
            hold all
            for sc = whichroi
                tubeROIs = CL.ROIgroups(sc).tubeROIs;
                for tc = 1:length(tubeROIs)
                    plot(tubeROIs(tc).interppath(:,1),tubeROIs(tc).interppath(:,2),'c.-')
                end
            end
            hold off
        end
    end
end