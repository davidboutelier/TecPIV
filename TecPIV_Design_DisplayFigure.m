function TecPIV_Design_DisplayFigure(hDisplayFigure, myhandles)
% Change display setting background UI

% Apply button at bottom right corner
hApplyDisplaySettingButton = uicontrol(...
    'Callback', @hApplyDisplaySetting,...
    'Style','pushbutton',...
    'String','Apply',...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[250 10 60 30]);

% show background image
hIMGcheckbox = uicontrol(...
    'Style','checkbox', ...
    'String', 'show background image', ...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'value', 1, ...
    'Position',[10 620 430 30]);

htextBackgroundColorPalette = uicontrol(...
    'Style','text', ...
    'String','Color map :',...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[10 585 100 30]);

hpopupBackgroundCPTSelector = uicontrol(...
    'Style','popupmenu',...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'String',{'parula','jet','hsv','hot','cool','spring','summer','autumn','winter','gray','bone','copper','pink','lines','colorcube','prism','flag','white'},...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[70 589 80 30]);

htextMinBackgroundColorPalette = uicontrol(...
    'Style','text', ...
    'String','Min :',...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[170 585 150 30]);

htextMaxBackgroundColorPalette = uicontrol(...
    'Style','text', ...
    'FontSize',8,...
    'String','Max :',...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[245 585 100 30]);

hMinBackgroundColorPalette = uicontrol(...
    'Style','edit',...
    'String','0',...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[195 596 35 20]);

hMaxBackgroundColorPalette = uicontrol(...
    'Style','edit',...
    'String','65536',...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[275 596 35 20]);

% show vectors
hveccheckbox = uicontrol(...
    'Style','checkbox', ...
    'String', 'show vectors', ...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'value', 1, ...
    'Position',[10 550 480 30]);

hbgvector = uibuttongroup(...
    'Units', 'pixels', ...
    'Parent', hDisplayFigure, ...
    'BorderType','None',...
    'Position',[10 305 150 250]);

hDisplayVectorRadioButton = uicontrol(...
    'Parent', hbgvector,...
    'Style','radiobutton',...
    'Value', 1, ...
    'String','display velocity',...
    'HorizontalAlignment','left',...
    'Units', 'pixels', ...
    'Position',[0 220 150 30]);

hDisplayGridRadioButton = uicontrol(...
    'Parent', hbgvector,...
    'Style','radiobutton',...
    'HorizontalAlignment','left',...
    'String','display grid',...
    'Units', 'pixels', ...
    'Position',[0 135 150 30]);

hDisplayStrainRadioButton = uicontrol(...
    'Parent', hbgvector,...
    'Style','radiobutton',...
    'HorizontalAlignment','left',...
    'String','display principal strain',...
    'Units', 'pixels', ...
    'Position',[0 75 250 30]);

htextVecScalingMode = uicontrol(...
    'Style','text', ...
    'String','Scaling mode:',...
    'HorizontalAlignment','left',...
    'FontSize',8,...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[10 490 150 30]);

hpopupVecColorScalingModeSelector = uicontrol(...
    'Style','popupmenu',...
    'String',{'mean','max', 'manual'},...
    'HorizontalAlignment','left',...
    'FontSize',8,...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[90 494 60 30]);

htextVecGridFactor = uicontrol(...
    'Style','text', ...
    'String','Grid resolution :',...
    'HorizontalAlignment','left',...
    'FontSize',8,...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[10 465 150 30]);

hpopupVecGridFactorSelector = uicontrol(...
    'Style','popupmenu',...
    'String',{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'},...
    'HorizontalAlignment','left',...
    'FontSize',8,...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[90 469 60 30]);

htextVecScalingLength = uicontrol(...
    'Style','text', ...
    'String','Reference value:',...
    'HorizontalAlignment','left',...
    'FontSize',8,...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[170 490 150 30]);

hVecScalingLength = uicontrol(...
    'Style','edit',...
    'String','0.1',...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[260 501 50 20]);

htextVecColor = uicontrol(...
    'Style','text', ...
    'String','Color:',...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[170 465 100 30]);

hpopupVecColorSelector = uicontrol(...
    'Style','popupmenu',...
    'String',{'yellow','magenta', 'cyan', 'red', 'green', 'blue', 'white', 'black'},...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[210 469 100 30]);

htextLagGridFactor = uicontrol(...
    'Style','text', ...
    'String','Grid resolution :',...
    'HorizontalAlignment','left',...
    'FontSize',8,...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[10 405 150 30]);

hpopupLagGridFactorSelector = uicontrol(...
    'Style','popupmenu',...
    'String',{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'},...
    'HorizontalAlignment','left',...
    'FontSize',8,...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[90 409 60 30]);

htextLagColor = uicontrol(...
    'Style','text', ...
    'String','Color:',...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[170 405 100 30]);

hpopupLagColorSelector = uicontrol(...
    'Style','popupmenu',...
    'String',{'yellow','magenta', 'cyan', 'red', 'green', 'blue', 'white', 'black'},...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[210 409 100 30]);

htexttypestrain = uicontrol(...
    'Style','text', ...
    'String','Type:',...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[10 345 100 30]);

hpopuptypestrainSelector = uicontrol(...
    'Style','popupmenu',...
    'String',{'extension', 'stretch'},...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[45 349 80 30]);

hshearcheckbox = uicontrol(...
    'Style','checkbox', ...
    'String', 'plot max shear', ...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'value', 0, ...
    'Position',[170 351 480 30]);

htextStrainScalingMode = uicontrol(...
    'Style','text', ...
    'String','Scaling mode:',...
    'HorizontalAlignment','left',...
    'FontSize',8,...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[10 320 150 30]);

hpopupStrainScalingModeSelector = uicontrol(...
    'Style','popupmenu',...
    'String',{'mean','max', 'manual'},...
    'HorizontalAlignment','left',...
    'FontSize',8,...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[90 324 60 30]);

htextStrainScalingLength = uicontrol(...
    'Style','text', ...
    'String','Reference value:',...
    'HorizontalAlignment','left',...
    'FontSize',8,...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[170 320 150 30]);

hStrainScalingLength = uicontrol(...
    'Style','edit',...
    'String','0.1',...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[260 331 50 20]);

htextStrainGridFactor = uicontrol(...
    'Style','text', ...
    'String','Grid resolution :',...
    'HorizontalAlignment','left',...
    'FontSize',8,...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[10 295 150 30]);

hpopupStrainGridFactorSelector = uicontrol(...
    'Style','popupmenu',...
    'String',{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'},...
    'HorizontalAlignment','left',...
    'FontSize',8,...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[90 299 60 30]);

htextStrainColor = uicontrol(...
    'Style','text', ...
    'String','Color:',...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[170 295 100 30]);

hpopupStrainColorSelector = uicontrol(...
    'Style','popupmenu',...
    'String',{'red-blue', 'blue-red','green-red', 'red-green'},...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[210 299 100 30]);

hscalarcheckbox = uicontrol(...
    'Style','checkbox', ...
    'String', 'show scalar', ...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'value', 1, ...
    'Position',[10 260 480 30]);

htextScalarType = uicontrol(...
    'Style','text', ...
    'String','Parameter:',...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'value', 1, ...
    'Position',[10 225 480 30]);

myhandles.ListScalarI = {'u', 'v', 'm', 'theta', 'du/dx', 'du/dy', 'dv/dx', 'dv/dy', 'exx', 'exy', 'eyy', 'omega', 'Ie', 'IIe', 'emax', 'emin', 'theta_p', 'gmax'};
myhandles.ListScalarF = {'Dx', 'Dy', 'm', 'theta', 'L', 'dDx/dx', 'dDx/dy', 'dDy/dx', 'dDy/dy', 'Exx', 'Exy', 'Eyy', 'Eyx', 'omega', 'IE', 'IIE', 'Emax', 'Emin', 'theta_p', 'Gmax'};

hpopupVecDerivativeTypeSelector = uicontrol(...
    'Style','popupmenu',...
    'String',myhandles.ListScalarI,...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'value', 1, ...
    'Position',[70 229 80 30]);

htextVecDerivativeDisplayRange = uicontrol(...
    'Style','text', ...
    'String','Range :',...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'value', 1, ...
    'Position',[10 201 480 30]);

hpopupVecDerivativeDisplayRangeSelector = uicontrol(...
    'Style','popupmenu',...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'String',{'min-max','+/- max','arbitrary range'},...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'value', 1, ...
    'Position',[70 205 80 30]);

hInterpDerivativetext = uicontrol(...
    'Style','text',...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'String','Interpolation :',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[170 201 70 30]);

hpopupVecDerivativeMethodSelector = uicontrol(...
    'Style','popupmenu',...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'String',{'linear','cubic','spline'},...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[240 205 70 30]);






htextScalarColorPalette = uicontrol(...
    'Style','text', ...
    'String','Color map :',...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[10 175 100 30]);

% put available cpt palettes in a list
CPT_Folder=fullfile(myhandles.TecPivFolder,'toolbox','colormaps');
DirEntries = dir(fullfile(CPT_Folder,'*.mat'));
ListOfCPTnames = {'parula','jet','hsv','hot','cool','spring','summer','autumn','winter','gray','bone','copper','pink','lines','colorcube','prism','flag','white'};
for Index = 1:length(DirEntries)
    CPTFileName = DirEntries(Index).name;
    [CPTfolder, CPTname, CPTextension] = fileparts(CPTFileName);
    ListOfCPTnames = [ListOfCPTnames strtrim(char(CPTname))];   
end  

ListOfCPTnames=char(ListOfCPTnames);

hpopupScalarCPTSelector = uicontrol(...
    'Style','popupmenu',...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'String',{ListOfCPTnames},...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[70 179 80 30]);

htextMinScalarColorPalette = uicontrol(...
    'Style','text', ...
    'String','Min :',...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[170 175 150 30]);

htextMaxScalarColorPalette = uicontrol(...
    'Style','text', ...
    'FontSize',8,...
    'String','Max :',...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[245 175 100 30]);

hMinScalarColorPalette = uicontrol(...
    'Style','edit',...
    'String','0',...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[195 186 35 20]);

hMaxScalarColorPalette = uicontrol(...
    'Style','edit',...
    'String','1',...
    'FontSize',8,...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'Position',[275 186 35 20]);


hroicheckbox = uicontrol(...
    'Style','checkbox', ...
    'String', 'show ROI', ...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'value', 1, ...
    'Position',[10 140 480 30]);

hroicheckbox = uicontrol(...
    'Style','checkbox', ...
    'String', 'show mask', ...
    'HorizontalAlignment','left',...
    'Parent', hDisplayFigure, ...
    'Units', 'pixels', ...
    'value', 1, ...
    'Position',[10 100 480 30]);














