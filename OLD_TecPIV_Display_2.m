function DisplayIMG(Ax,F, Datasets, Value,RawCpt,VectorField,Derivative)
%This function makes the figures. The image is not passed on to the function, only the dataset and frame number 
% F: the frame number to be displayed
% Datasets: the cell array with all the ddatasets in the current project
% Value: the number of the dataset to be employed
% RawCpt: a cell array with the colormaps
% VectorField: a cell array with the attributs to be used for the vectors
% Derivative: a cell array with the attributs to be used for the derivatives


    % Commun parameters regardless of data type
    TecPivFolder =  RawCpt{1,7};
    PIVStagePath=Datasets{Value,2}; % where the project was created (in the computer)
    ProjectName=Datasets{Value,3}; % name of project
    
    Dt = Datasets{Value,7}; % time increment for the dataset
    ValueDataImg = Datasets{Value,11}; % number of the dataset containing the images
    ImgPath = Datasets{ValueDataImg,1}; % path of the images
    
    DatasetName = Datasets{Value,1};

    % get the datatype from the datasets
    DataType=Datasets{Value,15};
    
    % attributs of the backgroundd image of the model 
    RawColorMapName=RawCpt{1,1};
    ColorMapMin=str2double(RawCpt{1,2});
    ColorMapMax=str2double(RawCpt{1,3});
    RelativeMin=ColorMapMin/65536;
    RelativeMax=ColorMapMax/65536;
    
    % scaling of image
    DoScale=RawCpt{1,4}; % do scale boolean
    ImScale=RawCpt{1,5}; % ImScale=6.84;
    Units=RawCpt{1,6};  % Units='mm';
    
    % select and clear the ax to make new figure
    axes(Ax);
    cla(Ax); % clear the axes
    
    switch DataType
        case 'image'
            % Plot image
            imgFramepath=fullfile(PIVStagePath,ProjectName,ImgPath,['IMG_' num2str(F) '.tif']);
            I=imread(imgFramepath);
    
            % adjust image to requested range
            I = imadjust(I,[RelativeMin RelativeMax],[0 1]);
            
            switch RawColorMapName
                    case 'jet'
                        colormap(jet);
                        subimage(I,jet(65536));
                    case 'gray'
                        colormap(gray);
                        subimage(I,gray(65536));
                end
                % add labels for positions in pix (scale is with scale bar)
                xlabel('pix');
                ylabel('pix');
            
            
            
        case 'vector'
             % Plot image
             % there is no option yet to plot the vectors and derivative
             % without the background image - Could be added
             
            imgFramepath=fullfile(PIVStagePath,ProjectName,ImgPath,['IMG_' num2str(F) '.tif']);
            I=imread(imgFramepath);
    
            % adjust image to requested range
            I = imadjust(I,[RelativeMin RelativeMax],[0 1]);
            
            % check if there is a mask
            if Datasets{Value,16} == 1 % there is a mask
                MaskDatasetNumber = Datasets{Value,17}; % dataset number where mask is
                MaskPath = Datasets{MaskDatasetNumber,1}; % path of the mask
            end 
            
            % Check if we display the vector field
            DisplayVector=VectorField{1,1};  % boolean display vecors yes no
            PlotVectorAsGrid = VectorField{1,16}; % boolean display as grid yes no
            VectorUnit = VectorField{1,9};
            VectorDensity = VectorField{1,2};
            VectorDisplayMode = VectorField{1,3};
            VectorScale = VectorField{1,4};
            
            if strcmp(VectorDisplayMode,'max') == 0 && strcmp(VectorDisplayMode,'mean') == 0
                VectorDisplayMode=VectorScale;
            end
            
            % check if we display the derivative
            DisplayDerivative=Derivative{1,1}; % boolean yes no
            DerivativeName=Derivative{1,2}; % name of the derivative
            RangeType=Derivative{1,3}; % Type of the data rane: minmax, maxmax, manual
            DeriveCpt=Derivative{1,4}; % name of the colormap
            
            if DisplayDerivative == 1
                
                % get the vectors
                X = VectorField{1,5};
                Y = VectorField{1,6};
                U = VectorField{1,7};
                V = VectorField{1,8};
        
                % get the minmax of vector field
                xmin=min2(X);
                xmax=max2(X);
                ymin=min2(Y);
                ymax=max2(Y);
                xwidth=[xmin xmax];
                ywidth=[ymin ymax];

                DX=abs(X(1,1)-X(1,2));
                DY=abs(Y(1,1)-Y(2,1));
                
                switch DerivativeName
                    case 'dUdx'
                        [dUdx,~] = gradient(U,DX);
                        DerivField = dUdx/(Dt);

                    case 'dVdy'
                        [~,dVdy] = gradient(V,DY);
                        DerivField = dVdy/(Dt);

                    case 'dUdy'
                         [~,dUdy] = gradient(U,DY);
                        DerivField = dUdy/(Dt);

                    case 'dVdx'
                        [dVdx,~] = gradient(V,DX);
                        DerivField = dVdx/(Dt);

                    case 'exx'
                        [dUdx,~] = gradient(U,DX);
                        DerivField = dUdx/(Dt);

                    case 'exy'
                        [~,dUdy] = gradient(U,DY);
                        [dVdx,~] = gradient(V,DX);
                        DerivField = 0.5*(dVdx/Dt + dUdy/Dt);

                    case 'eyx'
                        [~,dUdy] = gradient(U);
                        [dVdx,~] = gradient(V);
                        DerivField = 0.5*(dVdx/(DX*Dt) + dUdy/(DY*Dt));

                    case 'eyy'
                        [~,dVdy] = gradient(V);
                        DerivField = dVdy/(DY*Dt);

                    case 'vorticity'
                        [~,dUdy] = gradient(U,DY);
                        [dVdx,~] = gradient(V,DX);
                        DerivField = (dVdx - dUdy)/Dt;

                    case 'U'
                        DerivField = U;

                    case 'V'
                        DerivField = V;

                    case 'Theta'
                        [DerivField,~] = cart2pol(U,V);

                    case 'Rho'
                        [~,DerivField] = cart2pol(U,V);
                        DerivField = DerivField / Dt;

                    case 'Ie'
                        [dUdx,~] = gradient(U);
                        [~,dVdy] = gradient(V);
                        DerivField = dUdx / (DX*Dt) + dVdy / (DY*Dt); 

                    case 'IIe'
                        [dUdx,dUdy] = gradient(U);
                        [dVdx,dVdy] = gradient(V);
                        exy = 0.5* (dVdx/(DX*Dt) + dUdy/(DY*Dt));
                        eyx = exy;
                        exx = dUdx / (DX*Dt);
                        eyy = dVdy / (DY*Dt);
                        DerivField = exx.*eyy - exy.*eyx;    

                end
                
                InterpDeriv = Derivative{1,8}; % interpolate derivative yes or no
                
                if InterpDeriv == 1 
                    InterDerivMethod = Derivative{1,9} % how linear, spline...
                    switch InterDerivMethod
                        case 1
                            InterpMethod = 'linear';
                        case 2
                            InterpMethod = 'cubic';
                        case 3
                            InterpMethod = 'spline';
                    end
                    
                    
                    [HRX,HRY]=meshgrid(xmin:xmax,(ymin:ymax));
                    DerivField=interp2(X,Y,DerivField,HRX,HRY,InterDerivMethod);
                end
                
                % check that mask exist
                MaskExist = Datasets{Value,16};
                if MaskExist == 1 % if mask does exist
                    load(fullfile(PIVStagePath,ProjectName,MaskPath,['Mask_',num2str(F),'.mat']));
                    XM=RoiMask(:,1);
                    YM=RoiMask(:,2);
                    mask = poly2mask(XM-xmin,YM-ymin,size(DerivField,1),size(DerivField,2));
            
                    % delete extrapolated outside the mask area
                    DerivField(mask == 0) = 0;
                end
                
                ListMATLABCPT={'parula','jet','hsv','hot','cool','spring','summer',...
                'autumn','winter','gray','bone','copper','pink','lines',...
                'colorcube','prism','flag','white'};
                Lia = ismember(DeriveCpt,ListMATLABCPT); % check if colormap is in the list of MATLAB colormaps
                
                if Lia == 1 % it is a MATLAB colormap
                    RGB=DeriveCpt;
                    colormap(gca,RGB);
                else % it is a cutom colormap loaded from the folder
                    load(fullfile(TecPivFolder,'toolbox','colormaps',DeriveCpt));
                    colormap(gca,RGB);
                end 
                
                if RangeType == 1 % minmax            
                    MaxRange=max(max(DerivField));
                    MinRange=min(min(DerivField));
                    Range=[MinRange, MaxRange];

                elseif RangeType == 2 %+- max
                    AbsMaxRange=abs(max(max(DerivField)));
                    AbsMinRange=abs(min(min(DerivField)));
                    NewMax=max([AbsMaxRange, AbsMinRange]);
                    Range=[-NewMax,NewMax]

                else % manual mode
                    MinRange =Derivative{1,5};
                    MaxRange = Derivative{1,6};
                    Range = [MinRange, MaxRange];

                end
         
                DerivROI=imref2d(size(DerivField),xwidth,ywidth);
         
                % fuse background image with derivative
                I2=im2double(I); %convert image from uitnt16 to double
                IRGB = double2rgb(I2, gray); % convert to RGB
                
                DerivField(isnan(DerivField)) = 0; % remove NaN 
                size(DerivField)
                Ider=mat2im(DerivField,RGB,Range); % convert matrix to RGB

                RA=imref2d(size(I)); % area of image one (background)
                [D,~] = imfuse(IRGB,RA,Ider,DerivROI,'method','blend'); % blend the two images
                
                 subimage(D);
                iDeriv2=colorbar('location','westoutside');
                set(get(iDeriv2,'child'),'YData',Range);
                set(iDeriv2,'YLim',Range);  
                caxis(Range);
                set(gca, 'TickDir', 'out')
                hold on
                
                % check if we plot the ROI and the mask
                PlotROI = 1; %PlotROI = VectorField{1,X};
                PlotMask = 1; %PlotMask = VectorField{1,X};
                
                if PlotROI == 1
                    ROIX=[ROI(1), ROI(1)+ROI(3), ROI(1)+ROI(3), ROI(1),ROI(1)];
                    ROIY=[ROI(2), ROI(2), ROI(2)+ROI(4), ROI(2)+ROI(4),ROI(2)];
                    plot(ROIX,ROIY, 'r', 'LineWidth', 1)
                    hold on
                end
        
                if PlotMask == 1
                    plot(RoiMask(:,1),RoiMask(:,2), 'g', 'LineWidth', 1)
                    hold on
                end
                
            else  % no we do not plot the derivative of the vector field
                switch RawColorMapName
                    case 'jet'
                        colormap(jet);
                        subimage(I,jet(65536));
                    case 'gray'
                        colormap(gray);
                        subimage(I,gray(65536));
                end
                % add labels for positions in pix (scale is with scale bar)
                xlabel('pix');
                ylabel('pix');
                
   
            end
            
            if DisplayVector == 1  
                if PlotVectorAsGrid == 1
                else
                    % get the vectors
                    X = VectorField{1,5};
                    Y = VectorField{1,6};
                    U = VectorField{1,7};
                    V = VectorField{1,8};
                    typevector = VectorField{1,15};
                    
                    % if vector unit is pix, then no scaling and velocity
                    % is in pix/frame. No imscale, no dt
                    if strcmp(VectorUnit,'pix') == 0 
                        ImScale=VectorField{1,12};
                        U=U/(ImScale*Dt);
                        V=V/(ImScale*Dt);
                    end
                    
                    % make copy of UV
                    U1 = U; V1=V;
                    
                    % where vector is not good (typevector not 1), set UV to NaN
                    U1(typevector ~= 1) = NaN;
                    V1(typevector ~= 1) = NaN;
                    
                    % make copy of UV
                    U3= U; V3 =V;
                    
                    % where data is not interpolated set to NaN. Keep only interpolated vectors
                    % in patch
                    U3(typevector == 1) = NaN;
                    V3(typevector == 1) = NaN;
                    
                    % downsample UV
                    nx=VectorDensity(1);
                    ny=VectorDensity(2);
                    
                    XS = X(1:nx:end, 1:ny:end);
                    YS = Y(1:nx:end, 1:ny:end);
            
                    US3 = U3(1:nx:end, 1:ny:end);
                    VS3 = V3(1:nx:end, 1:ny:end);
            
                    US1 = U1(1:nx:end, 1:ny:end);
                    VS1 = V1(1:nx:end, 1:ny:end);

                    XS=XS(1,:);
                    YS=(YS(:,1))';
            
                    TF = ~any(~isnan(US3(:)));
            
                    % plot interpolated vectors
                    if TF == 0 % U3 is not empty
                        ncquiverref(XS,YS,US3,VS3,VectorUnit,VectorDisplayMode,'',[1 1 1],2); 
                        hold on
                    end
                    
                    % plot not interpolated vectors
                    ncquiverref(XS,YS,US1,VS1,VectorUnit,VectorDisplayMode,'true',[0.0 0.0 0.0],2);
                    hold on

                end
            end
        case 'Lagrangian'
            % load images
            imgFramepath=fullfile(PIVStagePath,ProjectName,ImgPath,['IMG_' num2str(F) '.tif']);
            I=imread(imgFramepath);
    
            % adjust image to requested range
            I = imadjust(I,[RelativeMin RelativeMax],[0 1]);
            
            % laod the lagrangian sum vectors
            load(fullfile(PIVStagePath,ProjectName, DatasetName, ['Vector_',num2str(F),'.mat']));
            
            %load the masks
            MaskDatasetNumber = Datasets{Value,17}; % dataset number where mask is
            MaskPath = Datasets{MaskDatasetNumber,1};
            load(fullfile(PIVStagePath,ProjectName,MaskPath,['Mask_',num2str(F),'.mat']));
            
            % patch the holes in grids
            XX=inpaint_nans(X,3);
            YY=inpaint_nans(Y,3);
            UU=inpaint_nans(U,3);
            VV=inpaint_nans(V,3);
            
            % get the minmax of vector field
            c = 0;
            if c == 1
                xmin=0;
                xmax=size(I,2)-1;
                ymin=0;
                ymax=size(I,1)-1;
            else
                xmin = floor(min2(XX));
                ymin = floor(min2(YY));
                xmax = ceil(max2(XX));
                ymax = ceil(max2(YY));
            end
            
            xwidth=[xmin xmax];
            ywidth=[ymin ymax];
            
            DX=abs(XX(1,1)-XX(1,2));
            DY=abs(YY(1,1)-YY(2,1));
            
            % define new x and y vectors
            Nx=xmin:1:xmax;
            Ny=ymin:1:ymax;

            % create new x and y grids
            [NX,NY]=meshgrid(Nx,Ny);
            
            % calculate gradient of lagrangian displacements
            % points are initially every 64 pixels (PIV resolution for each experiment)
            [exx,exy] = gradient(UU,DX,DY);
            [eyx,eyy] = gradient(VV,DX,DY);
            
            % interpolate the gradients to get finer resolution
            Nexx = griddata(XX,YY,exx,NX,NY,'linear');
            Nexy = griddata(XX,YY,exy,NX,NY,'linear');
            Neyy = griddata(XX,YY,eyy,NX,NY,'linear');
            Neyx = griddata(XX,YY,eyx,NX,NY,'linear');
            
            Derive = 0.5*(exy+eyx +(exx.*exy+eyx.*eyy));
            NDerive = 0.5*(Nexy+Neyx+(Nexx.*Nexy+Neyx.*Neyy));
            
            
            % check that mask exist
            MaskExist = Datasets{Value,16};
            
            if MaskExist == 1 % if mask does exist
                load(fullfile(PIVStagePath,ProjectName,MaskPath,['Mask_',num2str(F),'.mat']));
                XM=RoiMask(:,1);
                YM=RoiMask(:,2);
                mask = poly2mask(XM-xmin,YM-ymin,size(NDerive,1),size(NDerive,2));
                
                % delete extrapolated outside the mask area
                NDerive(mask == 0) = 0;
            end
            

            MaxRange=max(max(Derive));
            MinRange=min(min(Derive));
            Range=[MinRange, MaxRange];
            
            DerivROI=imref2d(size(NDerive),xwidth,ywidth);
            
            nx=4; ny=4;
            XS = X(1:nx:end, 1:ny:end);
            YS = Y(1:nx:end, 1:ny:end);
            US = U(1:nx:end, 1:ny:end);
            VS = V(1:nx:end, 1:ny:end);
            %typevectorS = typevector(1:nx:end, 1:ny:end);

            XS=XS(1,:);
            YS=(YS(:,1))';
            
            I2=im2double(I); %convert image from uitnt16 to double
            IRGB = double2rgb(I2, gray); % convert to RGB
            Ider=mat2im(NDerive,parula,Range); % convert matrix to RGB
            
            RA=imref2d(size(I)); % area of image one (background)
            [D,RD] = imfuse(IRGB,RA,Ider,DerivROI,'method','blend'); % blend the two images
            
            subimage(D);
            iDeriv2=colorbar('location','westoutside');
            set(get(iDeriv2,'child'),'YData',Range);
            set(iDeriv2,'YLim',Range);  
            caxis(Range);
            set(gca, 'TickDir', 'out')
            
            % calculate the outline of the interpolated grid
            k=1;
            XX2=XX(1:k:end,1:k:end);
            YY2=YY(1:k:end,1:k:end);
            
            N = size(XX2,1);
            M = size(XX2,2);
            
            for i=1:M-1
                for j=1:N
                    VX = [XX2(j,i), XX2(j,i+1)];
                    VY = [YY2(j,i), YY2(j,i+1)];
                    plot(VX,VY,'-','color','black')
                    hold on
                end
            end

            % then draw vertical lines
            for i=1:M
                for j=1:N-1
                    VX = [XX2(j,i), XX2(j+1,i)];
                    VY = [YY2(j,i), YY2(j+1,i)];
                    plot(VX,VY,'-','color','black')
                    hold on
                end
            end
            daspect([1 1 1])
            
    
    end
        
        % add scale bar if toggle = phys
        if strcmp(DoScale,'phys') == 1
            % Get the current axis limits
            xlim=get(gca,'xlim'); xp1=xlim(1); xp2=xlim(2);
            ylim=get(gca,'ylim'); yp1=ylim(1); yp2=ylim(2);

            % scalebar approximately 1/10 of image width
            LengthBarPix=(xp2-xp1)/10; % initial length pix
            LengthBarPhys=round2(LengthBarPix/ImScale,1); % closest rounded length phys
            
            if LengthBarPhys > 100
                LengthBarPhys=round2(LengthBarPhys,100);
            elseif LengthBarPhys > 10
                LengthBarPhys=round2(LengthBarPhys,10);
            end
            
            LengthBarPix=LengthBarPhys*ImScale; % actual length in pix of rounded length
            reftext=[num2str(LengthBarPhys),' ',Units,' '];

            % set padding around the scale bar
            padx=diff(xlim)/100; 
            pady=diff(ylim)/100;

            % Set x position of scale bar
            xend=xp2-padx;
            xstart=xend-LengthBarPix;
            
            % Plot reference text in lower right hand corner
            ht=text(xstart,yp1+pady,reftext,'Visible','off','Parent',gca,'FontSize',8.5,...
            'VerticalAlignment','Bottom','HorizontalAlignment','Right');
            textextent=get(ht,'Extent');

            % Draw patch over area of vector key 
            xl=textextent(1)-padx;
            xr=xp2;
            yt=yp2;
            yb=yp2-(textextent(2)+textextent(4)+pady);

            patch([xl; xl; xr; xr],[yb; yt; yt; yb],[2; 2; 2; 2],'w', 'LineWidth',0.5,'Parent',gca);

            % Redraw reference text on top of patch
            text(xstart,(yb+yt)/2,2.1,reftext,'Parent',gca,'FontSize',8.5,...
                 'VerticalAlignment','Middle','HorizontalAlignment','Right');
            hold on

            % Set y position of reference vector
            yend=yp2-(textextent(2)+textextent(4)/2);
            ystart=yend;

            lx = [xstart xend];
            ly = [ystart yend];
            lz = 3*ones(size(ly));

            % Plot the scalebar
            ScaleBar=line(lx,ly,lz,'LineWidth',1,'Color','black','Parent',gca);
            uistack(ScaleBar,'top')
            hold on
            
        end        
end
    
    


    
    
            

            


    
    
                 


