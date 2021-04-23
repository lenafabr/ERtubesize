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
            % max distance between tubule and sheet
            opt.maxtubedist = inf;
            
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
            closeenough = true(length(tubeROIs),1);
            for tc = 1:length(tubeROIs)
                tube = tubeROIs(tc);
                Ftubes(tc) = sum(sum(tube.dilmask.*double(CL.imgmem-CL.bgROI.avg)));
                Ltubes(tc) = tube.L;%+2*tube.dilum*CL.pxperum;
                
                tubepath = tube.interppath;
                sheetbound = sheetROI.bound;

                if (~isinf(opt.maxtubedist))
                    % check if tube is close enough to sheet
                    mindist = inf;
                    for pc = 1:size(sheetbound,1)
                        dists = sqrt(sum((tubepath - sheetbound(pc,:)).^2,2));
                        mindist = min(mindist,min(dists));
                    end
                    closeenough(tc) = (mindist< opt.maxtubedist*CL.pxperum);
                end
            end
            
            Ltubes = Ltubes(closeenough);
            Ftubes = Ftubes(closeenough);
            
            if (opt.mintubelen>0)
                goodind = find(Ltubes>opt.mintubelen);
            else
                goodind = 1:length(Ltubes);
            end
            
            
            
            totFtube = sum(Ftubes(goodind));
            totLtube = sum(Ltubes(goodind));
            
            Rtube = totFtube/Fsheet*Asheet/totLtube/pi/CL.pxperum;            
        end
        
        function showROImasks(CL,options)
            
            opt.whichroi = 1:length(CL.ROIgroups);
            opt.showlum = false;
            opt.label = false;

            if (exist('options','var'))
                opt = copyStruct(options,opt);
            end            
            
            %% show all sheet and tube masks and paths
            totmask = zeros(CL.ImgSize);
            for sc= opt.whichroi
                sheetROI = CL.ROIgroups(sc).sheetROI;
                tubeROIs = CL.ROIgroups(sc).tubeROIs;
                totmask = totmask + sheetROI.erodemask;
                paths = {};
                for tc = 1:length(tubeROIs)
                    totmask = totmask + tubeROIs(tc).dilmask;
                end
            end
            if (opt.showlum)
                imshowpair(CL.imglum,totmask)
            else
                imshowpair(CL.imgmem,totmask)
            end
            hold all
            
            if (opt.label)
                for sc = opt.whichroi
                    cent = mean(CL.ROIgroups(sc).sheetROI.bound);
                    text(cent(1),cent(2),sprintf('%d',sc),'Color','w')
                end
            end

            for sc = opt.whichroi
                tubeROIs = CL.ROIgroups(sc).tubeROIs;
                for tc = 1:length(tubeROIs)
                    plot(tubeROIs(tc).interppath(:,1),tubeROIs(tc).interppath(:,2),'c.-')
                end
            end
            hold off
        end
    end
end