function TecPIV_Display_LagGrid(I,Ax,RawCpt,VectorField,Derivative,Dt,MaskExist,ROI,RoiMask,RGB,RangeType, X,Y,U,V,DX,L)
%UNTITLED Summary of this function goes here

fprintf('Updating figure. Please wait...') % this is because the building of the figure may take a lot of time
tstart = tic;
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

% check if we display the derivative
DisplayDerivative=Derivative{1,1}; % boolean yes no
DerivativeName=Derivative{1,2}; % name of the derivative
RangeType=Derivative{1,3}; % Type of the data rane: minmax, maxmax, manual
DeriveCpt=Derivative{1,4}; % name of the colormap

% patch the holes in grids
XX=inpaint_nans(X,3);
YY=inpaint_nans(Y,3);

% get area of deformed grid
xmin = floor(min2(XX));
ymin = floor(min2(YY));
xmax = ceil(max2(XX));
ymax = ceil(max2(YY));
xwidth=[xmin xmax];
ywidth=[ymin ymax];

DY=DX;
            
% define new x and y vectors
Nx=xmin:1:xmax;
Ny=ymin:1:ymax;

% create new x and y grids
[NX,NY]=meshgrid(Nx,Ny);

tstartscalar = tic;

switch DerivativeName
    case 'U'
        UU=inpaint_nans(U,3);
        Derive = UU;
        NDerive = griddata(XX,YY,UU,NX,NY,'linear');
    case 'V'
        VV=inpaint_nans(V,3);
        Derive = VV;
        NDerive = griddata(XX,YY,VV,NX,NY,'linear');
        
    case 'M'
        UU=inpaint_nans(U,3);
        VV=inpaint_nans(V,3);
        Derive = sqrt(VV.^2 + UU.^2);
        NDerive = griddata(XX,YY,Derive,NX,NY,'linear');
        
    case 'theta'
        UU=inpaint_nans(U,3);
        VV=inpaint_nans(V,3);
        Derive = atan2(VV,UU)*180/pi;
        NDerive = griddata(XX,YY,Derive,NX,NY,'linear');
        
    case 'L'
        Derive  = inpaint_nans(L,3);
        NDerive = griddata(XX,YY,Derive,NX,NY,'linear');
        
    case 'dU/dx'
        UU=inpaint_nans(U,3);
        DX
        [exx,~] = gradient(UU,DX);
        Nexx = griddata(XX,YY,exx,NX,NY,'linear');
        Derive = exx;
        GTitle = 'dU/dx';
        NDerive = Nexx;
        
    case 'dU/dy'
        UU=inpaint_nans(U,3);
        [~,exy] = gradient(UU,DX,DY);
        Nexy = griddata(XX,YY,exy,NX,NY,'linear');
        Derive = exy;
        GTitle = 'dU/dy';
        NDerive = Nexy;
        
    case 'dV/dx'
        VV=inpaint_nans(V,3);
        [eyx,~] = gradient(VV,DX,DY);
        Neyx = griddata(XX,YY,eyx,NX,NY,'linear');
        Derive = eyx;
        GTitle = 'dV/dx';
        NDerive = Neyx;
        
    case 'dV/dy'
        VV=inpaint_nans(V,3);
        [~,eyy] = gradient(VV,DX,DY);
        Neyy = griddata(XX,YY,eyy,NX,NY,'linear');
        Derive = eyy;
        GTitle = 'dV/dy';
        NDerive = Neyy;
        
    case 'Exy'
        UU=inpaint_nans(U,3);
        VV=inpaint_nans(V,3);
        [exx,exy] = gradient(UU,DX,DY);
        [eyx,eyy] = gradient(VV,DX,DY);
        Nexx = griddata(XX,YY,exx,NX,NY,'linear');
        Nexy = griddata(XX,YY,exy,NX,NY,'linear');
        Neyy = griddata(XX,YY,eyy,NX,NY,'linear');
        Neyx = griddata(XX,YY,eyx,NX,NY,'linear');
        
        Derive = 0.5*(exy+eyx +(exx.*exy+eyx.*eyy));
        GTitle = 'Exy = 0.5*(exy+eyx +(exx.*exy+eyx.*eyy))';
        NDerive = 0.5*(Nexy+Neyx+(Nexx.*Nexy+Neyx.*Neyy));
        
    case 'Exx'
        UU=inpaint_nans(U,3);
        VV=inpaint_nans(V,3);
        [exx,exy] = gradient(UU,DX,DY);
        [eyx,eyy] = gradient(VV,DX,DY);
        Nexx = griddata(XX,YY,exx,NX,NY,'linear');
        Nexy = griddata(XX,YY,exy,NX,NY,'linear');
        Neyy = griddata(XX,YY,eyy,NX,NY,'linear');
        Neyx = griddata(XX,YY,eyx,NX,NY,'linear');
        
        Derive = 0.5*(exx+exx +(exx.*exx+eyx.*eyx));
        GTitle = 'Exx = 0.5*(exx+exx +(exx.*exx+eyx.*eyx));';
        NDerive = 0.5*(Nexx+Nexx +(Nexx.*Nexx+Neyx.*Neyx));
    
    case 'Eyy'
        UU=inpaint_nans(U,3);
        VV=inpaint_nans(V,3);
        [exx,exy] = gradient(UU,DX,DY);
        [eyx,eyy] = gradient(VV,DX,DY);
        Nexx = griddata(XX,YY,exx,NX,NY,'linear');
        Nexy = griddata(XX,YY,exy,NX,NY,'linear');
        Neyy = griddata(XX,YY,eyy,NX,NY,'linear');
        Neyx = griddata(XX,YY,eyx,NX,NY,'linear');
        
        Derive = 0.5*(eyy+eyy +(exy.*exy+eyy.*eyy));
        GTitle = 'Eyy = 0.5*(eyy+eyy +(exy.*exy+eyy.*eyy))';
        NDerive = 0.5*(Neyy+Neyy +(Nexy.*Nexy+Neyy.*Neyy));
        
    case 'Eyx'
        UU=inpaint_nans(U,3);
        VV=inpaint_nans(V,3);
        [exx,exy] = gradient(UU,DX,DY);
        [eyx,eyy] = gradient(VV,DX,DY);
        Nexx = griddata(XX,YY,exx,NX,NY,'linear');
        Nexy = griddata(XX,YY,exy,NX,NY,'linear');
        Neyy = griddata(XX,YY,eyy,NX,NY,'linear');
        Neyx = griddata(XX,YY,eyx,NX,NY,'linear');
        
        Derive = 0.5*(eyx+exy +(eyy.*eyx+exy.*exx));
        GTitle = 'Eyx = 0.5*(eyx+exy +(eyy.*eyx+exy.*exx))';
        NDerive = 0.5*(Neyx+Nexy +(Neyy.*Neyx+Nexy.*Nexx));
    
    case 'IE'
        GTitle = 'IE';
        
        UU=inpaint_nans(U,3);
        VV=inpaint_nans(V,3);
        [exx,exy] = gradient(UU,DX,DY);
        [eyx,eyy] = gradient(VV,DX,DY);
        Nexx = griddata(XX,YY,exx,NX,NY,'linear');
        Nexy = griddata(XX,YY,exy,NX,NY,'linear');
        Neyy = griddata(XX,YY,eyy,NX,NY,'linear');
        Neyx = griddata(XX,YY,eyx,NX,NY,'linear');
        
        Exx = 0.5*(exx+exx +(exx.*exx+eyx.*eyx));
        Eyy = 0.5*(eyy+eyy +(exy.*exy+eyy.*eyy));
        Derive = Exx+Eyy;
        
        NExx = 0.5*(Nexx+Nexx +(Nexx.*Nexx+Neyx.*Neyx));
        NEyy = 0.5*(Neyy+Neyy +(Nexy.*Nexy+Neyy.*Neyy));
        NDerive = NExx+NEyy;
        
    case 'IIE'
        GTitle = 'IIE';
        
        UU=inpaint_nans(U,3);
        VV=inpaint_nans(V,3);
        [exx,exy] = gradient(UU,DX,DY);
        [eyx,eyy] = gradient(VV,DX,DY);
        Nexx = griddata(XX,YY,exx,NX,NY,'linear');
        Nexy = griddata(XX,YY,exy,NX,NY,'linear');
        Neyy = griddata(XX,YY,eyy,NX,NY,'linear');
        Neyx = griddata(XX,YY,eyx,NX,NY,'linear');
        
        Exx = 0.5*(exx+exx +(exx.*exx+eyx.*eyx));
        Eyy = 0.5*(eyy+eyy +(exy.*exy+eyy.*eyy));
        Exy = 0.5*(exy+eyx +(exx.*exy+eyx.*eyy));
        Eyx = 0.5*(eyx+exy +(eyy.*eyx+exy.*exx));
        Derive = Exx.*Eyy - Exy.*Eyx;
        
        NExx = 0.5*(Nexx+Nexx +(Nexx.*Nexx+Neyx.*Neyx));
        NEyy = 0.5*(Neyy+Neyy +(Nexy.*Nexy+Neyy.*Neyy));
        NExy = 0.5*(Nexy+Neyx +(Nexx.*Nexy+Neyx.*Neyy));
        NEyx = 0.5*(Neyx+Nexy +(Neyy.*Neyx+Nexy.*Nexx));
        NDerive = NExx.*NEyy - NExy.*NEyx;
         
end
tendscalar = toc(tstartscalar);

tstartmask = tic;
% delete extrapolated outside the mask area
XM=RoiMask(:,1);
YM=RoiMask(:,2);
mask = poly2mask(XM - floor(min2(X(~isnan(X)))),YM -floor(min2(Y(~isnan(Y)))),size(NDerive,1),size(NDerive,2));
NDerive(mask == 0) = 0;
NDerive(isnan(NDerive) == 1) = 0;
tendsmask = toc(tstartmask);

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

GridCol = char(VectorField{1,14});
k=VectorField{1,2};
k= k(1);

if k == 1
    startx =1;
    starty =1;
else
    SX = size(XX);
    SSX = SX(2);
    SSY = SX(1);
    
    rX = mod(SSX-1,k);
    rY = mod(SSY-1,k);
    
    if rX == 0
        startx = 1;
    else
        if mod(rX,2) == 0
            startx = 1 + rX/2;
        else
            startx = 1 + (rX - 1)/2;
        end    
    end

    if rY == 0
        starty = 1;
    else
        if mod(rY,2) == 0
            starty = 1 + rY/2;
        else
            starty = 1 + (rY - 1)/2;
        end   
    end    
end

XX2=XX(starty:k:end,startx:k:end);
YY2=YY(starty:k:end,startx:k:end);

% new custom solution
N = size(XX2,1);
M = size(XX2,2);

colormap(gca,RGB);
I = double2rgb(im2double(I),gray);
Ider = mat2im(DerivField,RGB,Range);
DerivROI=imref2d(size(DerivField),xwidth,ywidth);
    
RA = imref2d(size(I));

tstartfuse = tic;
[D,~] = imfuse(I,RA,Ider,DerivROI,'method','blend');
tendfuse = toc(tstartfuse);
    
subimage(D)
iDeriv2=colorbar('location','westoutside');
set(get(iDeriv2,'child'),'YData',Range);
set(iDeriv2,'YLim',Range);  
caxis(Range);
set(gca, 'TickDir', 'out')
ylabel(iDeriv2,DerivativeName);
hold on 


for i = 1:N
    line(XX2(i,1:end),YY2(i,1:end),'color',GridCol);
    hold on
end

for j = 1:M
    line(XX2(1:end,j),YY2(1:end,j),'color',GridCol);
end


% for i=1:M-1
%     for j=1:N
%         VX = [XX2(j,i), XX2(j,i+1)];
%         VY = [YY2(j,i), YY2(j,i+1)];
%         plot(VX,VY,'-','color',GridCol)
%         hold on
%     end
% end
% 
% % then draw vertical lines
% for i=1:M
%     for j=1:N-1
%         VX = [XX2(j,i), XX2(j+1,i)];
%         VY = [YY2(j,i), YY2(j+1,i)];
%         plot(VX,VY,'-','color',GridCol)
%         hold on
%     end
% end

daspect([1 1 1])

tstartdrawmask = tic;
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

tenddrawmask = toc(tstartdrawmask);

%% If image is scaled. Add a scalebar
tstartdrawscale = tic;
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
tendscale = toc(tstartdrawscale);

tend = toc(tstart);
fprintf(' done.\n')
fprintf('- calculating the scalar field in %f s\n',tendscalar)
fprintf('- mask filter scalar plot in %f s\n', tendsmask)
fprintf('- bitmaps fused in %f s\n', tendfuse)
fprintf('- drawing ROI and mask in %f s\n', tenddrawmask)
fprintf('- drawing scale bar in %f s\n', tendscale)
fprintf('=> figure created in %f s\n', tend)

end

