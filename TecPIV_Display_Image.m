function TecPIV_Display_Image(I,Ax,RawCpt,VectorField,Derivative)
%TecPIV_Display_Image produce the figure when the display is image type

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

% display image
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

