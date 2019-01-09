function TecPIV_Display_Vectors(I,Ax,RawCpt,VectorField,Derivative,Dt,MaskExist,ROI,RoiMask,RGB,RangeType)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% select and clear the ax to make new figure
axes(Ax);
cla(Ax); % clear the axes

% attributs of the backgroundd image of the model 
RawColorMapName=RawCpt{1,1};  % name of colormap
ColorMapMin=str2double(RawCpt{1,2});  % value for min
ColorMapMax=str2double(RawCpt{1,3});  % value for max
RelativeMin=ColorMapMin/65536;
RelativeMax=ColorMapMax/65536;

% scaling of image
DoScale=RawCpt{1,4}; % do scale boolean
ImScale=RawCpt{1,5}; % ImScale eg 6.84 pix/mm;
Units=RawCpt{1,6};  % Units eg 'mm';

% adjust image to requested range
I = imadjust(I,[RelativeMin RelativeMax],[0 1]);

% Check if we display the vector field
DisplayVector=VectorField{1,1};  % boolean display vecors yes no
PlotVectorAsGrid = VectorField{1,16}; % boolean display as grid yes no
VectorUnit = VectorField{1,9};
VectorDensity = VectorField{1,2};
VectorDisplayMode = VectorField{1,3};
VectorScale = VectorField{1,4};

VectorCol = char(VectorField{1,14});

if strcmp(VectorDisplayMode,'max') == 0 && strcmp(VectorDisplayMode,'mean') == 0
    VectorDisplayMode=VectorScale;
end

% check if we display the derivative
DisplayDerivative=Derivative{1,1}; % boolean yes no
DerivativeName=Derivative{1,2}; % name of the derivative
RangeType=Derivative{1,3}; % Type of the data rane: minmax, maxmax, manual
DeriveCpt=Derivative{1,4}; % name of the colormap

%% make fused background image if plotting derivative 
if DisplayDerivative == 1
    % get the vectors
    X = VectorField{1,5};
    Y = VectorField{1,6};
    U = VectorField{1,7};
    V = VectorField{1,8};
    
    InterpMethod = 3;
    XX = inpaint_nans(X,InterpMethod);
    YY = inpaint_nans(Y,InterpMethod);
    UU = inpaint_nans(U,InterpMethod);
    VV = inpaint_nans(V,InterpMethod); 
        
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
        case 'du/dx'
            [dUdx,~] = gradient(UU,DX);
            DerivField = dUdx/(Dt);

        case 'dv/dy'
            [~,dVdy] = gradient(VV,DY);
            DerivField = dVdy/(Dt);

        case 'du/dy'
            [~,dUdy] = gradient(UU,DY);
            DerivField = dUdy/(Dt);

        case 'dv/dx'
            [dVdx,~] = gradient(VV,DX);
            DerivField = dVdx/(Dt);

        case 'exx'
            [dUdx,~] = gradient(UU,DX);
            DerivField = dUdx/(Dt);

        case 'exy'
            [~,dUdy] = gradient(UU,DY);
            [dVdx,~] = gradient(VU,DX);
            DerivField = 0.5*(dVdx/Dt + dUdy/Dt);

        case 'eyx'
            [~,dUdy] = gradient(UU);
            [dVdx,~] = gradient(VV);
            DerivField = 0.5*(dVdx/(DX*Dt) + dUdy/(DY*Dt));

        case 'eyy'
            [~,dVdy] = gradient(VV);
            DerivField = dVdy/(DY*Dt);

        case 'omega'
            [~,dUdy] = gradient(UU,DY);
            [dVdx,~] = gradient(VV,DX);
            DerivField = (dVdx - dUdy)/Dt;

        case 'u'
            DerivField = UU;

        case 'v'
            DerivField = VV;

        case 'theta'
            [DerivField,~] = cart2pol(UU,VV);

        case 'm'
            [~,DerivField] = cart2pol(UU,VV);
            DerivField = DerivField / Dt;

        case 'Ie'
            [dUdx,~] = gradient(UU);
            [~,dVdy] = gradient(VV);
            DerivField = dUdx / (DX*Dt) + dVdy / (DY*Dt); 

        case 'IIe'
            [dUdx,dUdy] = gradient(UU);
            [dVdx,dVdy] = gradient(VV);
            exy = 0.5* (dVdx/(DX*Dt) + dUdy/(DY*Dt));
            eyx = exy;
            exx = dUdx / (DX*Dt);
            eyy = dVdy / (DY*Dt);
            DerivField = exx.*eyy - exy.*eyx; 
            
        case 'theta_p'
            [dUdx,dUdy] = gradient(UU,DX,DY);
            [dVdx,dVdy] = gradient(VV,DX,DY);
            exx = dUdx/Dt;
            eyy = dVdy/Dt;
            exy = 0.5*(dVdx/Dt + dUdy/Dt);
            
            thetap = 0.5*atan(2*exy ./ (exx-eyy));
            DerivField = thetap;
            
        case 'emin'
            [dUdx,dUdy] = gradient(UU,DX,DY);
            [dVdx,dVdy] = gradient(VV,DX,DY);
            exx = dUdx/Dt;
            eyy = dVdy/Dt;
            exy = 0.5*(dVdx/Dt + dUdy/Dt);
            emin = 0.5*(exx+eyy) - ((0.5*(exx-eyy)).^2 + exy.^2).^0.5;
            DerivField = emin;
        
        case 'emax'
            [dUdx,dUdy] = gradient(UU,DX,DY);
            [dVdx,dVdy] = gradient(VV,DX,DY);
            exx = dUdx/Dt;
            eyy = dVdy/Dt;
            exy = 0.5*(dVdx/Dt + dUdy/Dt);
            emax = 0.5*(exx+eyy) + ((0.5*(exx-eyy)).^2 + exy.^2).^0.5;
            DerivField = emax;
            
        case 'gmax'
            [dUdx,dUdy] = gradient(UU,DX,DY);
            [dVdx,dVdy] = gradient(VV,DX,DY);
            exx = dUdx/Dt;
            eyy = dVdy/Dt;
            exy = 0.5*(dVdx/Dt + dUdy/Dt);
            emax = 0.5*(exx+eyy) + ((0.5*(exx-eyy)).^2 + exy.^2).^0.5;
            emin = 0.5*(exx+eyy) - ((0.5*(exx-eyy)).^2 + exy.^2).^0.5;
            gmax = emax - emin;
            DerivField = gmax;      
    end
    
    InterDerivMethod = Derivative{1,9}; % how linear, spline...
        switch InterDerivMethod
            case 1
                InterpMethod = 'linear';
            case 2
                InterpMethod = 'cubic';
            case 3
                InterpMethod = 'spline';
        end
     
    [HRX,HRY]=meshgrid(xmin:xmax,(ymin:ymax));
    NDerive=griddata(XX,YY,DerivField,HRX,HRY,InterDerivMethod);
    
    XM=RoiMask(:,1);
    YM=RoiMask(:,2);
    mask = poly2mask(XM - floor(min2(X(~isnan(X)))),YM - floor(min2(Y(~isnan(Y)))),size(NDerive,1),size(NDerive,2));
    NDerive(mask == 0) = 0;
    NDerive(isnan(NDerive) == 1) = 0;
    DerivField = NDerive;
    
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
    
    colormap(gca,RGB);
    I = double2rgb(im2double(I),gray);
    Ider = mat2im(DerivField,RGB,Range);
    DerivROI=imref2d(size(DerivField),xwidth,ywidth);
    
    RA = imref2d(size(I));
    [D,~] = imfuse(I,RA,Ider,DerivROI,'method','blend');
    
    subimage(D)
    iDeriv2=colorbar('location','westoutside');
    set(get(iDeriv2,'child'),'YData',Range);
    set(iDeriv2,'YLim',Range);  
    caxis(Range);
    set(gca, 'TickDir', 'out')
    ylabel(iDeriv2,DerivativeName);
    hold on 
    
else % otherwise use the backround image of model surface
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

% check if we plot the ROI and the mask
if MaskExist == 1
    
    PlotROI = VectorField{1,17};
    PlotMask = VectorField{1,18};
                
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
end

if DisplayVector == 1
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
    
    % downsample UV
	nx=VectorDensity(1);
	ny=VectorDensity(2);
                    
	XS = X(1:nx:end, 1:ny:end);
	YS = Y(1:nx:end, 1:ny:end);
            
	US = U(1:nx:end, 1:ny:end);
	VS = V(1:nx:end, 1:ny:end);
        
	XS=XS(1,:);
	YS=(YS(:,1))';
    
	% plot not interpolated vectors
	ncquiverref(XS,YS,US,VS,VectorUnit,VectorDisplayMode,'true',VectorCol,2);
	hold on
end

%% If image is scaled. Add a scalebar
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

