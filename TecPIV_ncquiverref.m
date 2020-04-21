function TecPIV_ncquiverref(x,y,u,v,units,reftype,refvec,veccol,cont)
% modified from ncquiverref written by Andrew Roberts

ScaleFontSize = 24;
ColorbarFontSize = 20;
VectorLineWidth = 2;

% make sure the mapping toolbox is present
h=ver('map'); 
if isempty(h) 
    error('Mapping toolbox not installed');
end

% don't color the vectors by length
col=false;

% get current axis 
h=get(gcf,'CurrentAxes');

% use meshgrid if needed
sx=size(x); sy=size(y);
if (sx(1)==1 && sy(1)==1) || (sx(2)==1 && sy(2)==1)
 [x,y]=meshgrid(x,y);
 
elseif sx(1)==1 || sx(2)==1
 error('Dimensions of x and y are inconsistent')
 
elseif sy(1)==1 || sy(2)==1
 error('Dimensions of x and y are inconsistent')
 
elseif sx~=sy
 error('Dimensions of x and y are inconsistent')
 
end

% check sizes all agree in input data
if size(x)~=size(y)
    error('x and y sizes disagree')
end

if size(x)~=size(u)
    error('x and u sizes disagree')
end

if size(y)~=size(v)
    error('y and v sizes disagree')
end


% get the magnitude of the vector field
[th,z] = cart2pol(u,v);

% remove masked grid points from the input by filling coordinates with NaN;
x(isnan(u))=NaN;
y(isnan(u))=NaN;

% Scale the vectors according to the reference arrow vector length based on
% the mean distance between grid points. This is a good measure, as it remains 
% constant for multiple plots using the same grid with different values.
x1=abs(diff(x')); x2=abs(diff(x)); 
y1=abs(diff(y')); y2=abs(diff(y));
[~,z1] = cart2pol(x1,y1); [~,z2] = cart2pol(x2,y2);
scalelength=0.8*min(mean(z1(~isnan(z1))),mean(z2(~isnan(z2))));

% Calculate reference vector length based on rounded median
% or maximum value of plot.  The default is median based.
if isnumeric(reftype) && ~col
	%disp('Calculating reference vector based on input number');
	refval=reftype;
    
elseif strcmp(reftype,'median') && ~col
	%disp('Calculating reference vector based on median');
    z(z==0)=NaN;
	refval=median(z(~isnan(z)))
    
elseif strcmp(reftype,'max') && ~col
	%disp('Calculating reference vector based on maximum');
	refval=max(z(~isnan(z)));
    
elseif ~col
	error('reftype must be either "max" or "median"');
end

% Remove NaN values that will not be plotted
% and turn points into a row of coordinates
u=u(~isnan(x))';
v=v(~isnan(x))';
y=y(~isnan(x))';
x=x(~isnan(x))';

% Set arrow size (1= full length of vector)
arrow=0.40;

% set scale value based on refval and scale length
roundp=floor(log10(refval));

refval=floor(refval/(10.^roundp))*(10.^roundp);
scale=scalelength/refval;

% Center vectors over grid points
xstart=x-0.5*scale*u;
xend=x+0.5*scale*u;
ystart=y-0.5*scale*v;
yend=y+0.5*scale*v;

% Get x coordinates of each vector plotted
lx = [xstart; x; ...
      xstart+(1-arrow/3)*(xend-xstart); ...
      xend-arrow*(scale*u+arrow*(scale*v)); ...
      xend; ...
      xend-arrow*(scale*u-arrow*(scale*v)); ...
      xstart+(1-arrow/3)*(xend-xstart); ...
      repmat(NaN,size(x))];

% Get y coordinates of each vector plotted
ly = [ystart; y; ...
      ystart+(1-arrow/3)*(yend-ystart); ...
      yend-arrow*(scale*v-arrow*(scale*u)); ...
      yend; ...
      yend-arrow*(scale*v+arrow*(scale*u)); ...
      ystart+(1-arrow/3)*(yend-ystart); ...
      repmat(NaN,size(y))];

% Plot the vectors
line(lx,ly,'Color',veccol, 'LineWidth', VectorLineWidth);

% Draw the reference vector key at altitude 2 above the map and grid
if refvec

    % Get the reference text string, formatted to powers of ten if required
    if refval < 0.1 | refval > 100 
        factor=floor(log10(refval));
        reftext=[num2str(refval/(10^factor)),' \times 10^{',num2str(factor),'} ',units,' '];
    else
        reftext=[num2str(refval),' ',units,' '];
    end
    
    % Get the current axis limits
    xlim=get(gca,'xlim'); xp1=xlim(1); xp2=xlim(2);
    ylim=get(gca,'ylim'); yp1=ylim(1); yp2=ylim(2);
    
    % set padding around the reference vector
    %padx=diff(xlim)/100; 
    %pady=diff(ylim)/100;
    
    padx = 50;
    pady = 10;
    
    % Set x position of reference vector
    xend=xp2-padx;
    xstart=xend-scalelength;
    
    % Plot reference text in lower right hand corner
	ht=text(xstart,yp1+pady,reftext,'Visible','off','Parent',gca,'FontSize',ScaleFontSize,...
        'VerticalAlignment','Bottom','HorizontalAlignment','Right');
	textextent=get(ht,'Extent');

	% Draw patch over area of vector key 
	xl=textextent(1)+0.5*padx;
	xr=xp2;
	yb=yp1;
	yt=textextent(2)+textextent(4)+pady;
	hp=patch([xl; xl; xr; xr],[yb; yt; yt; yb],[2; 2; 2; 2],'w',...
          'LineWidth',0.5,'Parent',gca);
	uistack(hp,'top');
    
    % Redraw reference text on top of patch
	ht=text(xstart,(yb+yt)/2,2.1,reftext,'Parent',gca,'FontSize',ScaleFontSize,...
         'VerticalAlignment','Middle','HorizontalAlignment','Right');
     
    % Set y position of reference vector
	yend=textextent(2)+textextent(4)/2;
	ystart=yend;

	% Get x coordinates of reference vector plotted
	lx = [xstart; ...
      xstart+(1-arrow/3)*(xend-xstart); ...
      xend-arrow*scalelength; ...
      xend; ...
      xend-arrow*scalelength; ...
      xstart+(1-arrow/3)*(xend-xstart); ...
      NaN];

	% Get y coordinates of reference vector plotted
	ly = [ystart; ...
      ystart+(1-arrow/3)*(yend-ystart); ...
      yend+arrow*(arrow*scalelength); ...
      yend; ...
      yend-arrow*(arrow*scalelength); ...
      ystart+(1-arrow/3)*(yend-ystart); ...
      NaN];

	% Get z coordinates of reference vector
	lz = 2*ones(size(ly));

	% Plot the reference vector
	line(lx,ly,lz,'Color',veccol,'LineWidth', VectorLineWidth);
    
end 
 

end