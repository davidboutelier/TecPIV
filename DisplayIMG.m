function DisplayIMG(F, Datasets, Value,RawCpt,VectorField,Derivative)
%This function makes the figures. The image is not passed on to the function, only the dataset and frame number 
% F: the frame number to be displayed
% Datasets: the cell array with all the ddatasets in the current project
% Value: the number of the dataset to be employed
% RawCpt: a cell array with the colormaps
% VectorField: a cell array with the attributs to be used for the vectors
% Derivative: a cell array with the attributs to be used for the derivatives

    % get the datatype from the datasets
    DataType=Datasets{Value,15};
    PIVStagePath=Datasets{Value,2}; % where the project was created (in the computer)
    ProjectName=Datasets{Value,3}; % name of project
    MaskDatasetNumber = Datasets{Value,17};
    MaskPath = Datasets{MaskDatasetNumber,1};
    Dt = Datasets{Value,7};
    ValueDataImg = Datasets{Value,11}; % number of the dataset containing the images
    ImgPath = Datasets{ValueDataImg,1};
    
    axes(Ax);
    cla(Ax); % clear the axes
    
    % attributs of the backgroundd image of the model 
    RawColorMapName=RawCpt{1,1};
    ColorMapMin=str2double(RawCpt{1,2});
    ColorMapMax=str2double(RawCpt{1,3});
    RelativeMin=ColorMapMin/65536;
    RelativeMax=ColorMapMax/65536;

    DoScale=RawCpt{1,4};
    ImScale=RawCpt{1,5}; %ImScale=6.84;
    Units=RawCpt{1,6}; %Units='mm';
    
    % Check if we display the vector field
    DisplayVector=VectorField{1,1};  % boolean display vecors yes no
    PlotVectorAsGrid = VectorField{1,16}; % boolean display as grid yes no

    DisplayDerivative=Derivative{1,1};
    DerivativeName=Derivative{1,2};
    RangeType=Derivative{1,3};
    DeriveCpt=Derivative{1,4};

    %DerivAlpha = Derivative{1,7};
    InterpDeriv = Derivative{1,8}; % interpolate derivative yes or no
    InterDerivMethod = Derivative{1,9}; % how linear, spline...
    
    ListMATLABCPT={'parula','jet','hsv','hot','cool','spring','summer',...
        'autumn','winter','gray','bone','copper','pink','lines',...
        'colorcube','prism','flag','white'};
    Lia = ismember(DeriveCpt,ListMATLABCPT);
    if Lia == 1
        RGB=DeriveCpt;
        colormap(gca,RGB);
    else
        RGB = load(fullfile(TecPivFolder,'toolbox','colormaps',DeriveCpt),'RGB');
        colormap(gca,RGB);
    end
    
    % if datasource is an image there is no vector to plot.
    if strcmp(DataType, 'image') == 1 
        DisplayVector = 0;
        PlotVectorAsGrid = 0;
    end

    % if data is incremental vector do not plot as grid
    if strcmp(DataType, 'vector') == 1
        PlotVectorAsGrid = 0;
    end
    
    %% Plot image
    imgFramepath=fullfile(PIVStagePath,ProjectName,ImgPath,['IMG_' num2str(F) '.tif']);
    I=imread(imgFramepath);
    
    % adjust image to requested range
    I = imadjust(I,[RelativeMin RelativeMax],[0 1]);
    
    % check if plot derive
    if DisplayDerivative == 1 % yes we do plot the derivative of vector field
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
        %xwidth=[xmin xmax];
        %ywidth=[ymin ymax];

        DX=abs(X(1,1)-X(1,2));
        DY=abs(Y(1,1)-Y(2,1));
        
        if InterpDeriv == 1
            [HRX,HRY]=meshgrid(xmin:xmax,(ymin:ymax));
        end
        
        switch DerivativeName
            case 'dUdx'
                [dUdx,~] = gradient(U);
                DerivField = dUdx/(DX*Dt);
                
            case 'dVdy'
                [~,dVdy] = gradient(V);
                DerivField = dVdy/(DY*Dt);
                
            case 'dUdy'
                 [~,dUdy] = gradient(U);
                DerivField = dUdy/(DY*Dt);
                
            case 'dVdx'
                [dVdx,~] = gradient(V);
                DerivField = dVdx/(DX*Dt);
                
            case 'exx'
                [dUdx,~] = gradient(U);
                DerivField = dUdx/(DX*Dt);
                
            case 'exy'
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
        
         if InterpDeriv == 1
            DerivField=interp2(X,Y,DerivField,HRX,HRY,InterDerivMethod);
         end
         
         % check that mask exist
        MaskExist = Datasets{Value,16};
        if MaskExist == 1 % if mask does exist
            
            RoiMask = load(fullfile(PIVStagePath,ProjectName,MaskPath,['Mask_',num2str(F),'.mat']),'RoiMask');
            
            XM=RoiMask(:,1);
            YM=RoiMask(:,2);
            mask = poly2mask(XM-xmin,YM-ymin,size(DerivField,1),size(DerivField,2));
        end
        
         % delete extrapolated outside the mask area
         DerivField(mask == 0) = 0;
         
         if RangeType == 1 % minmax            
             MaxRange=max(max(DerivField));
             MinRange=min(min(DerivField));
             Range=[MinRange, MaxRange];

         elseif RangeType == 2 %+- max
            AbsMaxRange=abs(max(max(DerivField)));
            AbsMinRange=abs(min(min(DerivField)));
            NewMax=max([AbsMaxRange, AbsMinRange]);
            Range=[-NewMax,NewMax];

         else % manual mode
           MinRange =Derivative{1,5};
           MaxRange = Derivative{1,6};
           Range = [MinRange, MaxRange];

         end
         
        
         % fuse background image with derivative
        I2=im2double(I); %convert image from uitnt16 to double
        IRGB = double2rgb(I2, gray); % convert to RGB
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
        
        if PlotROI == 1
            ROIX=[ROI(1), ROI(1)+ROI(3), ROI(1)+ROI(3), ROI(1),ROI(1)];
            ROIY=[ROI(2), ROI(2), ROI(2)+ROI(4), ROI(2)+ROI(4),ROI(2)];
            plot(ROIX,ROIY, 'g', 'LineWidth', 1)
            hold on
        end
        
        if PlotMask == 1
            plot(RoiMask(:,1),RoiMask(:,2), 'r', 'LineWidth', 1)
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
    
    if DisplayVector == 1 % we display the vectors
        % get the vectors
        X = VectorField{1,5};
        Y = VectorField{1,6};
        U = VectorField{1,7};
        V = VectorField{1,8};
        typevector = VectorField{1,15};

        if strcmp(VectorUnit,'pix') == 0 
            ImScale=VectorField{1,12};
            U=U/(ImScale*Dt);
            V=V/(ImScale*Dt);
        end
        
        if PlotVectorAsGrid % we display vectors as grid
        else % we display vectors as arrows
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
            
            % plot interpolated vectors
            ncquiverref(XS,YS,US3,VS3,'pix/frame',2,'',[0.25 0.25 0.25],2); 
            hold on
            % plot not interpolated vectors
            ncquiverref(XS,YS,US1,VS1,'pix/frame',2,'true',[0.0 0.0 0.0],2);
            hold on
            
            
        end
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
    end   
    
end

