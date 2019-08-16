% StandAloneFigureMaker

clear all
close all

%%% VARIABLES
F = 10;  % Frame
ImgPath = 'Y:\PIV\PROJECTS\MH5\Raw\Rectified';
VectorPath = 'Y:\PIV\PROJECTS\MH5\Raw\Rectified\Vectors';
ImScale = 6.73;
ColorMapMin = 0;
ColorMapMax = 20000;

DX = 16;
DY = 16;

Dt = 5;

Scalar = 'dudy';


cropXmin = 400;
cropXmax = 2500;
cropYmin = 1450;
cropYmax = 2650;

load(fullfile('C:\Users\dpb509\Documents\MATLAB\TecPIV','toolbox','colormaps','PurBluWhi.mat'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load tif model surface
RelativeMin=ColorMapMin/65536;
RelativeMax=ColorMapMax/65536;
ImRange = [RelativeMin RelativeMax];
Im = imread([ImgPath, '\IMG_', num2str(F),'.tif']);
Im2 = imadjust(Im,ImRange,[0 1]);

% load vectors
VectorData = load([VectorPath,'\Vector_',num2str(F),'.mat']);
X = VectorData.X;
Y = VectorData.Y;
U = VectorData.U;
V = VectorData.V;
Type = VectorData.typevector;
clear VectorData

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

switch Scalar
    case 'dudy'
        [~,dudy] = gradient(UU);
        DerivField = dudy / (DY*Dt);
        
    case 'IIe'
            [dUdx,dUdy] = gradient(UU);
            [dVdx,dVdy] = gradient(VV);
            exy = 0.5* (dVdx/(DX*Dt) + dUdy/(DY*Dt));
            eyx = exy;
            exx = dUdx / (DX*Dt);
            eyy = dVdy / (DY*Dt);
            DerivField = exx.*eyy - exy.*eyx; 
end

[HRX,HRY]=meshgrid(xmin:xmax,(ymin:ymax));
NDerive=griddata(XX,YY,DerivField,HRX,HRY,'linear');
%VectorData = load([VectorPath,'\Vector_',num2str(F),'.mat']);
RoiData = load([VectorPath,'\Mask_',num2str(F),'.mat']);
RoiMask=RoiData.RoiMask;

XM=RoiMask(:,1);
YM=RoiMask(:,2);
mask = poly2mask(XM - floor(min2(X(~isnan(X)))),YM - floor(min2(Y(~isnan(Y)))),size(NDerive,1),size(NDerive,2));
NDerive(mask == 0) = 0;
NDerive(isnan(NDerive) == 1) = 0;



% MaxRange= max2(DerivField);
MinRange= min2(DerivField);
%MinRange = -1e-6;
% NMax = max([abs(MinRange), abs(MaxRange)]);
% 
% NMax = 1e-4;
% Range=[-NMax, NMax];

Range = [MinRange 0];

DerivField = NDerive;

colormap(gca,RGB);
I = double2rgb(im2double(Im2),gray);
Ider = mat2im(DerivField,RGB,Range);
DerivROI=imref2d(size(DerivField),xwidth,ywidth);
    
RA = imref2d(size(I));
[D,~] = imfuse(I,RA,Ider,DerivROI,'method','blend');



U=U/(ImScale*Dt);
V=V/(ImScale*Dt);

% downsample UV
nx=4;
ny=4;
                    
XS = X(1:nx:end, 1:ny:end);
YS = Y(1:nx:end, 1:ny:end);          
US = U(1:nx:end, 1:ny:end);
VS = V(1:nx:end, 1:ny:end);       
XS=XS(1,:);
YS=(YS(:,1))';



subimage(D)
colormap(gca,RGB);
iDeriv2=colorbar('location','westoutside');
set(get(iDeriv2,'child'),'YData',Range);
set(iDeriv2,'YLim',Range);
%set(iDeriv2, 'Ticks', [-1e-4 -0.5e-4 0 0.5e-4 1e-4]);
caxis(Range);
set(gca, 'TickDir', 'out')
ylabel(iDeriv2,Scalar);
hold on 
axis([cropXmin cropXmax cropYmin cropYmax])
hold on
ncquiverref(XS,YS,US,VS,'mm/s',0.1,'true','black',2);

set(gca, 'Unit', 'pixels')

LbarPhy=50;
LbarPix = LbarPhy*ImScale

pad = 50;

xend = cropXmax - pad;
xstart = xend - LbarPix;
ystart = cropYmax -1*pad;
yend = cropYmax -1*pad;

lx = [xstart xend];
ly = [ystart yend];
lz = 3*ones(size(ly));

reftext=[num2str(LbarPhy),' mm '];
hold on
ht=text(xstart,ystart,3, reftext,'Visible','on','Parent',gca,'FontSize',8.5,...
    'VerticalAlignment','Middle','HorizontalAlignment','Right','Color','black');


textextent=get(ht,'Extent');

hold on
textextent=get(ht,'Extent');
% Draw patch over area of vector key 
xl= xstart - (textextent(3));
xr=cropXmax;
yt=cropYmax;
yb=cropYmax-(textextent(4)/2 +pad);
hp=patch([xl; xl; xr; xr],[yb; yt; yt; yb],[1; 1; 1; 1],'w', 'LineWidth',0.5,'Parent',gca);
hold on

% Plot the scalebar
ScaleBar=line(lx,ly,lz,'LineWidth',2,'Color','black','Parent',gca);
hold on

% set(gca, 'FontSize',12);
% g=get(gca,'Children')
% g = g([1 3 2 4:end])
% set(gca,'children',g)
% g=get(gca,'Children')


h = gcf();
% Set up the paper size / position
set(h, 'PaperUnits', 'centimeters');
set(h, 'PaperSize', [21 29.7]);    % Set to final desired size here as well as 2 lines below
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition', [0 0 21 29.7]);
exportname = fullfile('Y:\PIV\PROJECTS\MH5\',['export-',num2str(F),'.pdf']);
print(h,'-painters','-dpdf',exportname)
close(h)


