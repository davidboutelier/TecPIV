function TecPIV()
% TecPIV Summary of this function goes here
%   Detailed explanation goes here


%% Executes just before GUI is made visible.
function TecPIV_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject,handles);

end

%% Outputs from this function are returned to the command line.
function varargout = TecPIV_OutputFcn(hObject, ~, handles) 
guidata(hObject, handles);
varargout{1} = handles.output;
end


%% Create and then hide the GUI as it is being constructed.
screensize = get( groot, 'Screensize' );
hMainFigure = figure(...
    'name','main', ...
    'Position',[0,0,0.84*screensize(3),0.88*screensize(4)], ...
    'Visible','off', ...
    'MenuBar','none', ...
    'Toolbar','none', ...
    'Resize','off',...
    'CloseRequestFcn', @CloseRequestMainFigure, ...
    'HandleVisibility','callback', ...
    'Color', get(groot,...
    'defaultuicontrolbackgroundcolor'));

function CloseRequestMainFigure(hMainFigure,eventdata,handles)
    if isfield(myhandles,'ProjectID')
        warning('off','MATLAB:Figure:FigureSavedToMATFile')
        disp('Saving figures...')
        saveas(hMainFigure, 'MainFigure.m', 'm');
        saveas(hImportFigure, 'ImportFigure.m', 'm');
        
        
        saveas(hSecondFigure, 'SecondFigure.m', 'm');
        
        % create structure of handles
        Data = guidata(hMainFigure); 
               
        save(fullfile(myhandles.PathData, myhandles.ProjectID, myhandles.ProjectIDtag),'Data');
        cd(myhandles.TecPivFolder);
        
        s=isempty(gcp('nocreate'));
        if s == 0
            poolobj = gcp('nocreate');
            delete(poolobj);
        end
        
        disp('Closing now... bye bye')
        diary off
        pause(1);
        fclose('all');
        pause(1);
        eval(['diary ',logname]);
        copyfile(logname,fullfile(myhandles.PathData, myhandles.ProjectID, [ myhandles.ProjectID '.log'] ));
        fclose('all');
        pause(1);
        delete(findobj(0,'type','figure'));
    else
        delete(findobj(0,'type','figure'));
    end
end
function CloseRequestSecondFigure(hSecondFigure,eventdata,handles)
    set(get(hSecondFigure,'children'),'visible','off');
    hSecondFigure.Visible='off';
end


%% Create a import figure 
hImportFigure = figure(...
    'name','import', ...
    'Position',[0,0,210,100], ...
    'Visible','off', ...
    'MenuBar','none', ...
    'Toolbar','none', ...
    'HandleVisibility','callback', ...
    'Resize','off',...
    'CloseRequestFcn', @CloseRequestImportFigure, ...
    'Color', get(groot,...
    'defaultuicontrolbackgroundcolor'));

hCPFigure = figure(...
    'name', 'control points', ...
    'Position',[0,0,260,200], ...
    'Visible','off', ...
    'MenuBar','none', ...
    'Toolbar','none', ...
    'HandleVisibility','callback', ...
    'Resize','off',...
    'CloseRequestFcn', @CloseRequestImportFigure, ...
    'Color', get(groot,...
    'defaultuicontrolbackgroundcolor'));

%%
hSecondFigure = figure(...
    'name','TecPIV: settings', ...
    'Position',[0,0,0.5*screensize(3),0.66*screensize(4)], ...
    'Visible','off', ...
    'MenuBar','none', ...
    'Toolbar','none', ...
    'HandleVisibility','callback', ...
    'Resize','off',...
    'CloseRequestFcn', @CloseRequestSecondFigure, ...
    'Color', get(groot,...
    'defaultuicontrolbackgroundcolor'));

%% Create a display figure 
hDisplayFigure = figure(...
    'name','display', ...
    'Position',[0,0,340,650], ...
    'Visible','off', ...
    'MenuBar','none', ...
    'Toolbar','none', ...
    'HandleVisibility','callback', ...
    'Resize','off',...
    'CloseRequestFcn', @CloseRequestDisplayFigure, ...
    'Color', get(groot,...
    'defaultuicontrolbackgroundcolor'));



%% Create structure of handles
myhandles = guihandles(hMainFigure);


TecPivFile = which('TecPIV');
[TecPivFolder,name,ext] = fileparts(TecPivFile);
cd(TecPivFolder);

myhandles.TecPivFolder = TecPivFolder;

warning('off','MATLAB:dispatcher:nameConflict')
addpath(fullfile(TecPivFolder,'toolbox'));

% create a fontsize variable  - different for the platforms
if ismac()
    myhandles.MyFontSize = 10;
elseif isunix
    myhandles.MyFontSize = 8;       
else
    myhandles.MyFontSize = 9;
end

%% Create version
myhandles.Version='v.1801';
% date-version: 2.1.1 - 4/04/2016 --> v.1604
% date-version: 2.1.2 - 12/12/2016 -->v.1612
%% Create a log 
DATO=date;
logname=['matlab-log-',DATO];

disp(['This MATLAB session will be logged in file tecpiv-log-',DATO])
disp([' '])

warning('off','MATLAB:hg:uicontrol:StringMustBeNonEmpty')

ExistDiary=exist(logname,'file');

if ExistDiary == 2 %  means a file named 'diary' exist.
    delete(logname);
end

eval(['diary ',logname]);

% welcome message
disp(['Welcome. Starting a new session of TecPIV ' myhandles.Version ' ...'])

%% Open Parpool
s=isempty(gcp('nocreate'));

if s == 1 % matlabpool is not open
    
    % open new parallel pool with very long idle time
    %poolobj = parpool;
    %poolobj.IdleTimeout = 600000;
    
    parpool('IdleTimeout', 600000); 
    warning('off','MATLAB:datetime:InvalidSystemTimeZone')
    
end

%% Save the structure in hMainFigure
guidata(hMainFigure,myhandles); 

%% Query the compute capabilities of GPU devices.
fprintf('\n')
fprintf('Determine if computation can be performed on GPU....................')

f=fullfile(TecPivFolder ,'CUDA.mat');
load(f,'No_CUDA_msg');

if ismac
    % Code to run on Mac plaform
    fprintf(' No.\n')
    fprintf('\n')
    disp('CUDA is not yet made to work on macos.')
elseif isunix
    fprintf(' No.\n')
    fprintf('\n')
    disp('CUDA is not yet made to work on unix.')
else
    [status,cmdout] = system('nvcc');
    tf = strcmp(No_CUDA_msg,cmdout);
    if tf == 1 
        fprintf(' No.\n')
        fprintf('\n')
        disp('CUDA is not installed on this machine. Either there is no suitable GPU for CUDA compute, or CUDA is not installed properly.')

    else
    
        fprintf(' Yes.\n')
        fprintf('\n')
    
        g = gpuDevice;
        MaxComp=0;
    
    for ii = 1:gpuDeviceCount
        g = gpuDevice(ii);
        if g.ComputeCapability >= MaxComp
            MaxComp=g.ComputeCapability;
            Index=g.Index;
        end
        
    end
    
    gpuDevice(Index);
    fprintf('GPU Device %i has been selected\n', ...
       g.Index)
     
    end
    
end

%% //UI settings
% UI main figure
for folding=true
    % first we create the controls and the axe in the main window. 
    for folding2 = true
    
% get size MainWindow
P = hMainFigure.Position;

% Create player
hpanelPlayer = uipanel(...
    'Title','Controls',...
    'FontSize',myhandles.MyFontSize,...
    'Parent', hMainFigure, ...
    'Units', 'pixels', ...
    'Position',[10 10 200 (P(4) -20)]);
    
% Create the axes
hPlotAxes = axes(...    
    'Parent', hMainFigure, ...
    'Units', 'pixels', ...
    'HandleVisibility','callback', ...
    'Position',[220 10 (P(3)-220 -10) (P(4) -20) ]);

PP = hpanelPlayer.Position;

htextSourceSelector = uicontrol(...
    'Style','text', ...
    'FontSize', myhandles.MyFontSize,...
    'String','Select Data: ',...
    'Parent', hpanelPlayer, ...
    'Units', 'pixels', ...
    'HorizontalAlignment','left',...
    'Position',[5 (PP(4) - 60) 190 30]);

hpopupSourceSelector = uicontrol(...
    'Callback', @hSourceSelectorCallback,...
    'String',{'Cover'},...
    'FontSize',myhandles.MyFontSize,...
    'Style','popupmenu',...
    'Parent', hpanelPlayer, ...
    'Units', 'pixels', ...
    'HorizontalAlignment','left',...
    'Position',[5 (PP(4) - 80) 190 30]); 

htextPosition = uicontrol(...
    'Style','text', ...
    'FontSize', myhandles.MyFontSize,...
    'String','Position: ',...
    'Parent', hpanelPlayer, ...
    'Units', 'pixels', ...
    'HorizontalAlignment','left',...
    'Position',[5 (PP(4) - 120) 190 30]);

hslider = uicontrol(...
    'Callback', @hSliderCallback,...
    'Style','slider',...
    'Max',100,...
    'Min',1,...
    'Value',1,...
    'SliderStep',[0.05 0.2],...
    'Parent', hpanelPlayer, ...
    'Units', 'pixels', ...
    'Position',[5 (PP(4) - 130) 190 20]);

htextImgNum = uicontrol(...
    'Style','text', ...
    'FontSize',myhandles.MyFontSize,...
    'String','Number :', ...
    'Parent', hpanelPlayer, ...
    'Units', 'pixels', ...
    'HorizontalAlignment','left',...
    'Position',[5 (PP(4) - 175) 100 25]);

htextImgTime = uicontrol(...
    'Style','text', ...
    'FontSize',myhandles.MyFontSize,...
    'String','Time (s) :', ...
    'Parent', hpanelPlayer, ...
    'Units', 'pixels', ...
    'HorizontalAlignment','left',...
    'Position',[5 (PP(4) - 205) 100 25]);

hImgTime = uicontrol(...
    'Callback', @hImgTimeCallback,...
    'Style','edit',...
    'FontSize',myhandles.MyFontSize,...
    'String','0',...
    'Parent', hpanelPlayer, ...
    'Units', 'pixels', ...
    'HorizontalAlignment','left',...
    'Position',[70 (PP(4) - 171) 45 25]);

hImgNumber = uicontrol(...
    'Callback', @hImgNumCallback,...
    'Style','edit',...
    'FontSize',myhandles.MyFontSize,...
    'String','0',...
    'Parent', hpanelPlayer, ...
    'Units', 'pixels', ...
    'HorizontalAlignment','left',...
    'Position',[70 (PP(4) - 201) 45 25]);

    end
    
    % create the menus
    for folding3=true
    % Project menu
    hProjectMenu = uimenu(...       
        'Parent',hMainFigure,...
        'HandleVisibility','callback', ...
        'Label','Project');

    % Open Project menu item
    hOpenMenuitem = uimenu(...       
        'Parent',hProjectMenu,...
        'Label','Open',...
        'HandleVisibility','callback', ...
        'Callback', @hOpenMenuitemCallback);
    % New Project menu item
    hNewMenuitem = uimenu(...       
        'Parent',hProjectMenu,...
        'Label','New',...
        'HandleVisibility','callback', ...
        'Callback', @hNewMenuitemCallback);

    % Save & Close Project menu item
    hCloseMenuitem = uimenu(...       
        'Parent',hProjectMenu,...
        'Label','Save & Close',...
        'Separator','on',...
        'HandleVisibility','callback', ...
        'Callback', @hCloseMenuitemCallback);

    % Import menu
    hImportMenu = uimenu(...       
        'Parent',hMainFigure,...
        'HandleVisibility','callback', ...
        'Label','Import');

    % Import Experiment frames menu item
    hImportExpFramesMenuitem = uimenu(...       
        'Parent',hImportMenu,...
        'Label','Experiment frames',...
        'HandleVisibility','callback', ...
        'Callback', @hImportExpFramesMenuitemCallback);

    % Import Calibration frames menu item
    hImportCalibFramesMenuitem = uimenu(...       
        'Parent',hImportMenu,...
        'Label','Calibration frames',...
        'HandleVisibility','callback', ...
        'Callback', @hImportCalibFramesMenuitemCallback);

    % Calibrate menu
    hCalibrateMenu = uimenu(...       
        'Parent',hMainFigure,...
        'HandleVisibility','callback', ...
        'Label','Calibration');

    % Extract corners menu item
    hExtractControlPointsMenuitem = uimenu(...       
        'Parent',hCalibrateMenu,...
        'Label','Extract control points',...
        'HandleVisibility','callback', ...
        'Callback', @hExtractControlPointsMenuitemCallback);

    % urectify calibration menu item
    hRectifyCalibrationMenuitem = uimenu(...       
        'Parent',hCalibrateMenu,...
        'Label','Rectify calibration images',...
        'HandleVisibility','callback', ...
        'Callback', @hRectifyCalibrationMenuitemCallback);

    % Preprocess menu
    hPreprocessMenu = uimenu(...       
        'Parent',hMainFigure,...
        'HandleVisibility','callback', ...
        'Label','Preprocessing');

    hUndeformMenu = uimenu(...       
        'Parent',hPreprocessMenu,...
        'Label','Undeform dataset',...
        'HandleVisibility','callback', ...
        'Callback', @hUndeformMenuitemCallback);

    hContrastMenu = uimenu(...       
        'Parent',hPreprocessMenu,...
        'Label','Adjust contrast',...
        'HandleVisibility','callback', ...
        'Callback', @hContrastMenuCallback);

    % PIV processing menu
    hPIVprocessMenu = uimenu(...       
        'Parent',hMainFigure,...
        'HandleVisibility','callback', ...
        'Label','Processing');

    % PIV settings menu item
    hPIVsettingsMenuitem = uimenu(...       
        'Parent',hPIVprocessMenu,...
        'Label','PIV settings',...
        'HandleVisibility','callback', ...
        'Callback', @hPIVsettingsMenuitemCallback);

    % Display menu
    hDisplayMenu = uimenu(...       
        'Parent',hMainFigure,...
        'HandleVisibility','callback', ...
        'Label','Display');

    % display settings menu item
    hDisplaySettingsMenuitem = uimenu(...       
        'Parent',hDisplayMenu,...
        'Label','Settings',...
        'HandleVisibility','callback', ...
        'Callback', @hDisplaySettingsMenuitemCallback);

    % Export this frame
    hExportThisFrameMenuitem = uimenu(...       
        'Parent',hDisplayMenu,...
        'Label','Export this frame',...
        'HandleVisibility','callback', ...
        'Callback', @hExportThisFrameMenuitemCallback);

    % Export All frames
    hExportAllFrameMenuitem = uimenu(...       
        'Parent',hDisplayMenu,...
        'Label','Export all frames',...
        'HandleVisibility','callback', ...
        'Callback', @hExportAllFrameMenuitemCallback);

    % PIV postprocessing menu
    hPIVPostprocessMenu = uimenu(...       
        'Parent',hMainFigure,...
        'HandleVisibility','callback', ...
        'Label','Postprocessing');

    hLagrangianSumMenuitem = uimenu(...       
         'Parent',hPIVPostprocessMenu,...
         'Label','Calculate Lagrangian sum',...
         'HandleVisibility','callback', ...
         'Callback', @hLagrangianSumMenuitemCallback);

    % Eulerian sum
    hEulerianSumMenuitem = uimenu(...       
        'Parent',hPIVPostprocessMenu,...
        'Label','Calculate Eulerian sum',...
        'HandleVisibility','callback', ...
        'Callback', @hEulerianSumMenuitemCallback);
    end
end

% UI Display figure
for folding = true
    MyFontSize = myhandles.MyFontSize;
    
    P = hDisplayFigure.Position;
    
    % Apply button at bottom right corner
    hApplyDisplaySettingButton = uicontrol(...
        'Callback', @hApplyDisplaySetting,...
        'Style','pushbutton',...
        'String','Apply',...
        'FontSize',MyFontSize+2,...
        'HorizontalAlignment','left',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[(P(3) - 60) 10 60 30]);

    % show background image
    hIMGcheckbox = uicontrol(...
        'Style','checkbox', ...
        'String', 'show background image', ...
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'value', 1, ...
        'Position',[10 620 430 40]);

    htextBackgroundColorPalette = uicontrol(...
        'Style','text', ...
        'String','Color map :',...
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[10 585 100 30]);

    hpopupBackgroundCPTSelector = uicontrol(...
        'Style','popupmenu',...
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'String',{'parula','jet','hsv','hot','cool','spring','summer','autumn','winter','gray','bone','copper','pink','lines','colorcube','prism','flag','white'},...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[110 589 70 30]);

    htextMinBackgroundColorPalette = uicontrol(...
        'Style','text', ...
        'String','Min/Max:',...
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[200 585 150 30]);

    hMinBackgroundColorPalette = uicontrol(...
        'Style','edit',...
        'String','0',...
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[255 596 40 20]);

    hMaxBackgroundColorPalette = uicontrol(...
        'Style','edit',...
        'String','65536',...
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[300 596 40 20]);

    % show vectors
    hveccheckbox = uicontrol(...
        'Style','checkbox', ...
        'String', 'show vectors', ...
        'FontSize',MyFontSize,...
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
        'FontSize',MyFontSize,...
        'String','display velocity',...
        'HorizontalAlignment','left',...
        'Units', 'pixels', ...
        'Position',[0 220 150 30]);

    hDisplayGridRadioButton = uicontrol(...
        'Parent', hbgvector,...
        'FontSize',MyFontSize,...
        'Style','radiobutton',...
        'HorizontalAlignment','left',...
        'String','display grid',...
        'Units', 'pixels', ...
        'Position',[0 135 150 30]);

    hDisplayStrainRadioButton = uicontrol(...
        'Parent', hbgvector,...
        'FontSize',MyFontSize,...
        'Style','radiobutton',...
        'HorizontalAlignment','left',...
        'String','display principal strain',...
        'Units', 'pixels', ...
        'Position',[0 75 250 30]);

    htextVecScalingMode = uicontrol(...
        'Style','text', ...
        'String','Scaling mode:',...
        'HorizontalAlignment','left',...
        'FontSize',MyFontSize,...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[10 490 150 30]);

    hpopupVecColorScalingModeSelector = uicontrol(...
        'Style','popupmenu',...
        'String',{'mean','max', 'manual'},...
        'HorizontalAlignment','left',...
        'FontSize',MyFontSize,...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[110 494 70 30]);

    htextVecGridFactor = uicontrol(...
        'Style','text', ...
        'String','Grid sampling :',...
        'HorizontalAlignment','left',...
        'FontSize',MyFontSize,...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[10 465 150 30]);

    hpopupVecGridFactorSelector = uicontrol(...
        'Style','popupmenu',...
        'String',{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'},...
        'HorizontalAlignment','left',...
        'FontSize',MyFontSize,...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[110 469 70 30]);

    htextVecScalingLength = uicontrol(...
        'Style','text', ...
        'String','Reference:',...
        'HorizontalAlignment','left',...
        'FontSize',MyFontSize,...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[200 490 150 30]);

    hVecScalingLength = uicontrol(...
        'Style','edit',...
        'String','0.1',...
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[270 501 50 20]);

    htextVecColor = uicontrol(...
        'Style','text', ...
        'String','Color:',...
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[200 465 100 30]);

    hpopupVecColorSelector = uicontrol(...
        'Style','popupmenu',...
        'String',{'yellow','magenta', 'cyan', 'red', 'green', 'blue', 'white', 'black'},...
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[270 469 70 30]);

    htextLagGridFactor = uicontrol(...
        'Style','text', ...
        'String','Grid sampling :',...
        'HorizontalAlignment','left',...
        'FontSize',MyFontSize,...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[10 405 150 30]);

    hpopupLagGridFactorSelector = uicontrol(...
        'Style','popupmenu',...
        'String',{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'},...
        'HorizontalAlignment','left',...
        'FontSize',MyFontSize,...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[110 409 70 30]);

    htextLagColor = uicontrol(...
        'Style','text', ...
        'String','Color:',...
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[200 405 100 30]);

    hpopupLagColorSelector = uicontrol(...
        'Style','popupmenu',...
        'String',{'yellow','magenta', 'cyan', 'red', 'green', 'blue', 'white', 'black'},...
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[270 409 70 30]);

    htexttypestrain = uicontrol(...
        'Style','text', ...
        'String','Type:',...
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[10 345 100 30]);

    hpopuptypestrainSelector = uicontrol(...
        'Style','popupmenu',...
        'String',{'extension', 'stretch'},...
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[110 349 70 30]);

    hshearcheckbox = uicontrol(...
        'Style','checkbox', ...
        'FontSize',MyFontSize,...
        'String', 'plot max shear', ...
        'HorizontalAlignment','left',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'value', 0, ...
        'Position',[200 351 480 30]);

    htextStrainScalingMode = uicontrol(...
        'Style','text', ...
        'String','Scaling mode:',...
        'HorizontalAlignment','left',...
        'FontSize',MyFontSize,...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[10 320 150 30]);

    hpopupStrainScalingModeSelector = uicontrol(...
        'Style','popupmenu',...
        'String',{'mean','max', 'manual'},...
        'HorizontalAlignment','left',...
        'FontSize',MyFontSize,...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[110 324 70 30]);

    htextStrainScalingLength = uicontrol(...
        'Style','text', ...
        'String','Reference:',...
        'HorizontalAlignment','left',...
        'FontSize',MyFontSize,...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[200 320 150 30]);

    hStrainScalingLength = uicontrol(...
        'Style','edit',...
        'String','0.1',...
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[270 331 70 20]);

    htextStrainGridFactor = uicontrol(...
        'Style','text', ...
        'String','Grid sampling :',...
        'HorizontalAlignment','left',...
        'FontSize',MyFontSize,...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[10 295 150 30]);

    hpopupStrainGridFactorSelector = uicontrol(...
        'Style','popupmenu',...
        'String',{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'},...
        'HorizontalAlignment','left',...
        'FontSize',MyFontSize,...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[110 299 70 30]);

    htextStrainColor = uicontrol(...
        'Style','text', ...
        'String','Color:',...
        'FontSize',8,...
        'HorizontalAlignment','left',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[200 295 100 30]);

    hpopupStrainColorSelector = uicontrol(...
        'Style','popupmenu',...
        'String',{'red-blue', 'blue-red','green-red', 'red-green'},...
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[270 299 70 30]);

    hscalarcheckbox = uicontrol(...
        'Style','checkbox', ...
        'FontSize',MyFontSize,...
        'String', 'show scalar', ...
        'HorizontalAlignment','left',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'value', 1, ...
        'Position',[10 260 480 30]);

    htextScalarType = uicontrol(...
        'Style','text', ...
        'String','Parameter:',...
        'FontSize',MyFontSize,...
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
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'value', 1, ...
        'Position',[110 229 70 30]);

    htextVecDerivativeDisplayRange = uicontrol(...
        'Style','text', ...
        'String','Range :',...
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'value', 1, ...
        'Position',[10 201 480 30]);

    hpopupVecDerivativeDisplayRangeSelector = uicontrol(...
        'Style','popupmenu',...
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'String',{'min-max','+/- max','arbitrary range'},...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'value', 1, ...
        'Position',[110 205 70 30]);

    hInterpDerivativetext = uicontrol(...
        'Style','text',...
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'String','Interp:',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[200 201 70 30]);

    hpopupVecDerivativeMethodSelector = uicontrol(...
        'Style','popupmenu',...
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'String',{'linear','cubic','spline'},...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[270 205 70 30]);

    htextScalarColorPalette = uicontrol(...
        'Style','text', ...
        'String','Color map :',...
        'FontSize',MyFontSize,...
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
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'String',{ListOfCPTnames},...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[110 179 70 30]);

    htextMinScalarColorPalette = uicontrol(...
        'Style','text', ...
        'String','Min/Max:',...
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[200 175 200 30]);


    hMinScalarColorPalette = uicontrol(...
        'Style','edit',...
        'String','0',...
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[255 186 40 20]);

    hMaxScalarColorPalette = uicontrol(...
        'Style','edit',...
        'String','1',...
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'Position',[300 186 40 20]);


    hroicheckbox = uicontrol(...
        'Style','checkbox', ...
        'String', 'show ROI', ...
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'value', 1, ...
        'Position',[10 140 480 30]);

    hroicheckbox = uicontrol(...
        'Style','checkbox', ...
        'String', 'show mask', ...
        'FontSize',MyFontSize,...
        'HorizontalAlignment','left',...
        'Parent', hDisplayFigure, ...
        'Units', 'pixels', ...
        'value', 1, ...
        'Position',[10 100 480 30]);



end

% UI Import figure
for folding = true
    
	hImportFrames = uicontrol(...
        'Callback', @hImportFramesCallback,...
        'Style','pushbutton',...
        'String','Import',...
        'FontSize',myhandles.MyFontSize + 2,...
        'HorizontalAlignment','left',...
        'Parent', hImportFigure, ...
        'Units', 'pixels', ...
        'Position',[150 10 60 30]);
    
    htextTimeIncrement = uicontrol(...
        'Style','text', ...
        'HorizontalAlignment','left',...
        'String','Time increment (s):',...
        'FontSize',myhandles.MyFontSize,...
        'Parent', hImportFigure, ...
        'Units', 'pixels', ...
        'Position',[10 40 150 30]);
    
     htextImageFormat = uicontrol(...
        'Style','text', ...
        'String','Image format :',...
        'FontSize',myhandles.MyFontSize,...
        'HorizontalAlignment','left',...
        'Parent', hImportFigure, ...
        'Units', 'pixels', ...
        'Position',[10 80 150 30]);
    

    hTimeInc = uicontrol(...
        'Style','edit',...
        'String','1',...
        'FontSize',myhandles.MyFontSize,...
        'HorizontalAlignment','left',...
        'Parent', hImportFigure, ...
        'Units', 'pixels', ...
        'Position',[150 52 50 20]);
    

    hImgFormat = uibuttongroup(...
        'Title','',...
        'FontSize',myhandles.MyFontSize,...
        'Parent', hImportFigure, ...
        'Units', 'pixels', ...
        'Position',[110 80 200 50],...
        'BorderType','None');

    hradiobuttonImportDNG= uicontrol(...
        'Style','radiobutton',...
        'FontSize',myhandles.MyFontSize,...
        'HorizontalAlignment','left',...
        'String','DNG',...
        'Parent', hImgFormat, ...
        'Units', 'pixels', ...
        'Value',  1,...
        'Position',[0 0 50 50]);
    
    hradiobuttonImportTIF= uicontrol(...
        'Style','radiobutton',...
        'FontSize',myhandles.MyFontSize,...
        'HorizontalAlignment','left',...
        'String','TIF',...
        'Parent', hImgFormat, ...
        'Units', 'pixels', ...
        'Value',  0,...
        'Position',[65 0 50 50]);
% 
%     hradiobuttonImportExpTiff= uicontrol(...
%         'Style','radiobutton',...
%         'String','16-bit Tif',...
%         'Parent', hImgFormat, ...
%         'Units', 'normalized', ...
%         'Position',[0.33 0.0 0.66 1]);
% 
%     


    
end

% UI GCP Figure
for folding = true
    
    hStartControlPointsExtraction = uicontrol(...
        'Callback', @hControlPointsExtractionStartCallback,...
        'Style','pushbutton',...
        'String','Start',...
        'FontSize',myhandles.MyFontSize + 2,...
        'HorizontalAlignment','left',...
        'Parent', hCPFigure, ...
        'Units', 'pixels', ...
        'Position',[200 10 60 30]); 
    
    htextDefinePhysU = uicontrol(...
        'Style','text', ...
        'String','Physical unit :',...
        'FontSize',myhandles.MyFontSize,...
        'HorizontalAlignment','right',...
        'Parent', hCPFigure, ...
        'Units', 'pixels', ...
        'Position',[10 40 200 30]);

    htextDefineSizeSqY = uicontrol(...
        'Style','text', ...
        'FontSize',myhandles.MyFontSize,...
        'HorizontalAlignment','right',...
        'String','Length of squares along Y axis :',...
        'Parent', hCPFigure, ...
        'Units', 'pixels', ...
        'Position',[10 70 200 30]);
    
    htextDefineSizeSqX = uicontrol(...
        'Style','text', ...
        'FontSize',myhandles.MyFontSize,...
        'HorizontalAlignment','right',...
        'String','Length of squares along X axis :',...
        'Parent', hCPFigure, ...
        'Units', 'pixels', ...
        'Position',[10 100 200 30]);
    
    htextDefineNbSqY = uicontrol(...
        'Style','text', ...
        'FontSize',myhandles.MyFontSize,...
        'HorizontalAlignment','right',...
        'String','Number of squares along Y axis :',...
        'Parent', hCPFigure, ...
        'Units', 'pixels', ...
        'Position',[10 130 200 30]);
    
    htextDefineNbSqX = uicontrol(...
        'Style','text', ...
        'FontSize',myhandles.MyFontSize,...
        'HorizontalAlignment','right',...
        'String','Number of squares along X axis :',...
        'Parent', hCPFigure, ...
        'Units', 'pixels', ...
        'Position',[10 160 200 30]);
    
    hPhysU = uicontrol(...
        'Style','edit',...
        'FontSize',myhandles.MyFontSize,...
        'String','mm',...
        'Parent', hCPFigure, ...
        'Units', 'pixels', ...
        'Position',[220 52 40 20]);
    
    
    hSizeSqY = uicontrol(...
        'Style','edit',...
        'FontSize',myhandles.MyFontSize,...
        'String','15',...
        'Parent', hCPFigure, ...
        'Units', 'pixels', ...
        'Position',[220 82 40 20]);
    
    hSizeSqY = uicontrol(...
        'Style','edit',...
        'FontSize',myhandles.MyFontSize,...
        'String','15',...
        'Parent', hCPFigure, ...
        'Units', 'pixels', ...
        'Position',[220 112 40 20]);
    
    hNbSqY = uicontrol(...
        'Style','edit',...
        'FontSize',myhandles.MyFontSize,...
        'String','20',...
        'Parent', hCPFigure, ...
        'Units', 'pixels', ...
        'Position',[220 142 40 20]);
    
	hNbSqX = uicontrol(...
        'Style','edit',...
        'FontSize',myhandles.MyFontSize,...
        'String','30',...
        'Parent', hCPFigure, ...
        'Units', 'pixels', ...
        'Position',[220 172 40 20]);
    
end

% make figures visible/invisible
for folding = true
    
    % Move the GUI to the center of the screen.
    movegui(hMainFigure,'center');
    movegui(hSecondFigure,'center');
    movegui(hDisplayFigure,'center');
    movegui(hImportFigure,'center');
    movegui(hCPFigure,'center');
    
    
    % Make the GUI visible/invisible.
    hMainFigure.Visible = 'on';
    hSecondFigure.Visible = 'off';
    hDisplayFigure.Visible = 'off';
    hImportFigure.Visible = 'off';
    
    % make panels visible/invisible
    hpanelImportFrames.Visible = 'off';
    hpanelGetControlPoints.Visible = 'off';
    hpanelRectifySettings.Visible='off';
    hpanelUndeform.Visible='off';
    hpanelContrast.Visible='off'; 
    hpanelPIVSettings.Visible='off';
    hpanelExportAllSettings.Visible='off';
    
end

% define settings & defaults
for folding = true
    % Define RawCptProperties
    myhandles.RawCpt={'gray','0','65536','pix','1','pix',TecPivFolder}; % cpt name, min, max, toggle pix/phys, scale (1 if toggle=pix), unit label

    % Define Vector field display properties
    myhandles.VectorField{1,1} = 0; % display yes/no
    myhandles.VectorField{1,2} = [4 4]; % density vector field
    myhandles.VectorField{1,3} = 'max'; % display mode
    myhandles.VectorField{1,4} = 1; % refernce value when manual scaling
    myhandles.VectorField{1,5} = 0; 
    myhandles.VectorField{1,6} = 0; 
    myhandles.VectorField{1,7} = 0; 
    myhandles.VectorField{1,8} = 0; 
    myhandles.VectorField{1,9} = 'phys'; % unit
    myhandles.VectorField{1,10} = []; % time
    myhandles.VectorField{1,11} = []; % inc
    myhandles.VectorField{1,12} = []; % IMscale
    myhandles.VectorField{1,13} = []; % datapath
    myhandles.VectorField{1,14} = []; % color
    myhandles.VectorField{1,15} = []; % typevector
    myhandles.VectorField{1,16} = []; % plotasgrid
    
    % Define Derivative display properties
    myhandles.Derivative{1,1} = 0;
    myhandles.Derivative{1,2}= 'eyx';
    myhandles.Derivative{1,3} = 1; % minmax
    myhandles.Derivative{1,4} = 'DIV256';% default palette for derivative
    myhandles.Derivative{1,5} = 0; % default min
    myhandles.Derivative{1,6} = 1; % default max
    myhandles.Derivative{1,7} = 0.5; % alpha
    myhandles.Derivative{1,8} = 0; % Toggle interpolation 0 = no, 1 == yes
    myhandles.Derivative{1,9} = 'linear'; % default interpolation linear cubic or spline

    myhandles.ImScale=1; % in order to work with unscaled images
    myhandles.PhysU='pix'; % by default if no calibration unit is pix
    myhandles.NumberOfDatasets = 0; % no dataset at begining
    
    % Define defaults img
    hpopupBackgroundCPTSelector.Value= 10; % gray is default cpt for background image
    hMinBackgroundColorPalette.String='0';
    hMaxBackgroundColorPalette.String='65536';

    guidata(hMainFigure,myhandles);
    
end


%% //COVER
% place cover image in ax
filename = 'TecPIVCover.png';
y = imread(filename, 'BackgroundColor', [1 1 1]);
figure(hMainFigure);
Ax = hPlotAxes;
axes(Ax);
cla(Ax);
imshow(y, 'Parent', Ax);



for folding=true




% % Get control points
% hpanelGetControlPoints = uipanel(...
%     'Title','Get control points',...
%     'FontSize',8,...
%     'Parent', hSecondFigure, ...
%     'Units', 'normalized', ...
%     'Position',[0.01 0.01 0.98 0.98]);

% htextDefineCalibrationBoard = uicontrol(...
%     'Style','text', ...
%     'String','Calibration board for camera:',...
%     'FontWeight', 'bold', ...
%     'Parent', hpanelGetControlPoints, ...
%     'Units', 'normalized', ...
%     'HorizontalAlignment','left',...
%     'Position',[0.01 0.825 0.25 0.15]);

% hCamNumCalib = uicontrol(...
%     'Style','edit',...
%     'String','1',...
%     'Parent', hpanelGetControlPoints, ...
%     'Units', 'normalized', ...
%     'Position',[0.225 0.95 0.075 0.03]);


% 
% htextDefineNbSqX = uicontrol(...
%     'Style','text', ...
%     'String','Number of squares along X axis :',...
%     'Parent', hpanelGetControlPoints, ...
%     'Units', 'normalized', ...
%     'HorizontalAlignment','left',...
%     'Position',[0.01 0.775 0.225 0.15]);
% 
% htextDefineNbSqY = uicontrol(...
%     'Style','text', ...
%     'String','Number of squares along Y axis :',...
%     'Parent', hpanelGetControlPoints, ...
%     'Units', 'normalized', ...
%     'HorizontalAlignment','left',...
%     'Position',[0.01 0.725 0.225 0.15]);
% 
% htextDefineSizeSqX = uicontrol(...
%     'Style','text', ...
%     'String','Length of squares along X axis :',...
%     'Parent', hpanelGetControlPoints, ...
%     'Units', 'normalized', ...
%     'HorizontalAlignment','left',...
%     'Position',[0.35 0.775 0.225 0.15]);
% 
% htextDefineSizeSqY = uicontrol(...
%     'Style','text', ...
%     'String','Length of squares along Y axis :',...
%     'Parent', hpanelGetControlPoints, ...
%     'Units', 'normalized', ...
%     'HorizontalAlignment','left',...
%     'Position',[0.35 0.725 0.225 0.15]);

% hNbSqX = uicontrol(...
%     'Style','edit',...
%     'String','23',...
%     'Parent', hpanelGetControlPoints, ...
%     'Units', 'normalized', ...
%     'Position',[0.225 0.9 0.075 0.03]);
% 
% hNbSqY = uicontrol(...
%     'Style','edit',...
%     'String','15',...
%     'Parent', hpanelGetControlPoints, ...
%     'Units', 'normalized', ...
%     'Position',[0.225 0.85 0.075 0.03]);

% hSizeSqX = uicontrol(...
%     'Style','edit',...
%     'String','15',...
%     'Parent', hpanelGetControlPoints, ...
%     'Units', 'normalized', ...
%     'Position',[0.56 0.9 0.075 0.03]);
% 
% hSizeSqY = uicontrol(...
%     'Style','edit',...
%     'String','15',...
%     'Parent', hpanelGetControlPoints, ...
%     'Units', 'normalized', ...
%     'Position',[0.56 0.85 0.075 0.03]);



hStartControlPointsExtraction = uicontrol(...
    'Callback', @hControlPointsExtractionStartCallback,...
    'Style','pushbutton',...
    'String','Start',...
    'Parent', hpanelGetControlPoints, ...
    'Units', 'normalized', ...
    'Position',[0.75 0.89 0.1 0.05]);

hDoneControlPointsExtraction = uicontrol(...
    'Callback', @hControlPointsExtractionDoneCallback,...
    'Style','pushbutton',...
    'String','Done',...
    'Parent', hpanelGetControlPoints, ...
    'Units', 'normalized', ...
    'Position',[0.75 0.82 0.1 0.05]);

htextInstructionControlPointsTitle = uicontrol(...
    'Style','text', ...
    'String','Instructions :',...
    'FontWeight', 'bold', ...
    'Parent', hpanelGetControlPoints, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','left',...
    'Position',[0.01 0.65 0.15 0.15]);

hInstructionControlPoints = uicontrol(...
    'Style','Text',...
    'Parent', hpanelGetControlPoints, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','left',...
    'Position',[0.01 0.70 0.98 0.6]);
InstructionMessage = {'Click on the four extreme corners on the rectangular checkerboard pattern.',...
          'The first point is the origin point of the reference frame attached to the grid. The other three points of the rectangular grid can be clicked in any order.'};
% Wrap string
[outstring,newpos]=textwrap(hInstructionControlPoints,InstructionMessage);
hInstructionControlPoints.String=outstring;
hInstructionControlPoints.Position=newpos;

% Rectification settings UI
hpanelRectifySettings = uipanel(...
    'Title','Rectification settings',...
    'FontSize',8,...
    'Parent', hSecondFigure, ...
    'Units', 'normalized', ...
    'Position',[0.01 0.01 0.98 0.98]);

htextFrameNum = uicontrol(...
    'Style','text', ...
    'String','Use calibration image number :',...
    'Parent', hpanelRectifySettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','left',...
    'Position',[0.01 0.775 0.225 0.15]);

hRectFrameNum = uicontrol(...
    'Style','edit',...
    'String','1',...
    'Parent', hpanelRectifySettings, ...
    'Units', 'normalized', ...
    'Position',[0.225 0.9 0.075 0.03]);

htextRectMethodSelector = uicontrol(...
    'Style','text', ...
    'String','Rectification method :',...
    'Parent', hpanelRectifySettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','left',...
    'Position',[0.01 0.7275 0.225 0.15]);

hpopupRectMethodSelector = uicontrol(...
    'Style','popupmenu',...
    'Parent', hpanelRectifySettings, ...
    'String',{'Projective','Polynomial','Projective + polynomial'},...
    'Units', 'normalized', ...
    'Position',[0.225 0.733 0.220 0.15]);

htextOrderPoly = uicontrol(...
    'Style','text', ...
    'String','Order of polynomial function :',...
    'Parent', hpanelRectifySettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','left',...
    'Position',[0.01 0.675 0.220 0.15]);

hOrderPoly = uicontrol(...
    'Style','edit',...
    'String','3',...
    'Parent', hpanelRectifySettings, ...
    'Units', 'normalized', ...
    'Position',[0.225 0.8 0.075 0.03]);

hStartRectifyCalib = uicontrol(...
    'Callback', @hStartRectifyCalibCallback,...
    'Style','pushbutton',...
    'String','Start',...
    'Parent', hpanelRectifySettings, ...
    'Units', 'normalized', ...
    'Position',[0.75 0.89 0.1 0.05]);

hDoneRectifyCalib = uicontrol(...
    'Callback', @hDoneRectifyCalibCallback,...
    'Style','pushbutton',...
    'String','Done',...
    'Parent', hpanelRectifySettings, ...
    'Units', 'normalized', ...
    'Position',[0.75 0.82 0.1 0.05]);

htextInstructionRectifyTitle = uicontrol(...
    'Style','text', ...
    'String','Instructions :',...
    'FontWeight', 'bold', ...
    'Parent', hpanelRectifySettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','left',...
    'Position',[0.01 0.60 0.15 0.15]);

hInstructionRectify = uicontrol(...
    'Style','Text',...
    'Parent', hpanelRectifySettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','left',...
    'Position',[0.01 0.65 0.98 0.6]);

InstructionMessage = {'Specifies the order of polynomials to use. order can be 2, 3, or 4. The higher the order of the polynomial, the better the fit, but the result can contain more curves than the base image. Minimum Number of Control Point Pairs: 6 (order 2), 10 (order 3) and 15 (order 4)'};
% Wrap string
[outstring,newpos]=textwrap(hInstructionRectify,InstructionMessage);
hInstructionRectify.String=outstring;
hInstructionRectify.Position=newpos;

% Undeform all images UI
hpanelUndeform = uipanel(...
    'Title','Undeform dataset',...
    'FontSize',8,...
    'Parent', hSecondFigure, ...
    'Units', 'normalized', ...
    'Position',[0.01 0.01 0.98 0.98]);

hradiobuttonRotateImages= uicontrol(...
    'Style','radiobutton',...
    'String','Rotate images',...
    'Parent', hpanelUndeform, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','left',...
    'Position',[0.01 0.85 0.3 0.05]);

htextDefineAngleRotation = uicontrol(...
      'Style','edit',...
    'String','0',...
    'Parent', hpanelUndeform, ...
    'Units', 'normalized', ...
    'Position',[0.175 0.855 0.075 0.03]);

hStartUndeform = uicontrol(...
    'Callback', @hStartUndeformCallback,...
    'Style','pushbutton',...
    'String','Start',...
    'Parent', hpanelUndeform, ...
    'Units', 'normalized', ...
    'Position',[0.75 0.89 0.1 0.05]);

% Adjust Contrast UI
hpanelContrast = uipanel(...
    'Title','Adjust contrast',...
    'FontSize',8,...
    'Parent', hSecondFigure, ...
    'Units', 'normalized', ...
    'Position',[0.01 0.01 0.98 0.98]);

hInverseImageRadioButton = uicontrol(...
    'Style','radiobutton',...
    'String','Inverse image',...
    'Parent', hpanelContrast, ...
    'Units', 'normalized', ...
    'Position',[0.01 0.95 0.2 0.025]);

hSubtractBackgroundRadioButton = uicontrol(...
    'Style','radiobutton',...
    'String','Subtract background',...
    'Parent', hpanelContrast, ...
    'Units', 'normalized', ...
    'Position',[0.01 0.9 0.2 0.025]);

hNormalizeIntensityRadioButton = uicontrol(...
    'Style','radiobutton',...
    'String','Normalize intensity',...
    'Parent', hpanelContrast, ...
    'Units', 'normalized', ...
    'Position',[0.01 0.85 0.2 0.025]);

hUseMaskRadioButton = uicontrol(...
    'Style','radiobutton',...
    'String','Use mask',...
    'Parent', hpanelContrast, ...
    'Units', 'normalized', ...
    'Position',[0.01 0.80 0.2 0.025]);

hGaussianWindowSizeText = uicontrol(...
    'Style','edit',...
    'String','50',...
    'Parent', hpanelContrast, ...
    'Units', 'normalized', ...
    'Position',[0.4 0.9 0.075 0.05]);

hGaussianSigmaText = uicontrol(...
    'Style','edit',...
    'String','10',...
    'Parent', hpanelContrast, ...
    'Units', 'normalized', ...
    'Position',[0.4 0.8 0.075 0.05]);

hApplyContrast = uicontrol(...
    'Callback', @hApplyContrastCallback,...
    'Style','pushbutton',...
    'String','Start',...
    'Parent', hpanelContrast, ...
    'Units', 'normalized', ...
    'Position',[0.6 0.7 0.1 0.05]);

% PIV settings UI
hpanelPIVSettings = uipanel(...
    'Title','PIV settings',...
    'FontSize',8,...
    'Parent', hSecondFigure, ...
    'Units', 'normalized', ...
    'Position',[0.01 0.01 0.98 0.98]);

htextDatasetPIV = uicontrol(...
    'Style','text', ...
    'String','Dataset:',...
    'FontWeight', 'bold', ...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','left',...
    'Position',[0.01 0.825 0.25 0.15]);

% hpopupSourceSelector.String
hpopupPIVDataSelector = uicontrol(...
    'Style','popupmenu',...
    'String',{'item 1', 'item2'},...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'Position',[0.075,0.835,0.35,0.15]);

htextPIVImgStart = uicontrol(...
    'Style','text', ...
    'String','Start image :',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','right',...
    'Position',[0.45 0.825 0.125 0.15]);

htextPIVImgEnd = uicontrol(...
    'Style','text', ...
    'String','End image :',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','right',...
    'Position',[0.45 0.775 0.125 0.15]);

htextPIVImgInc = uicontrol(...
    'Style','text', ...
    'String','Increment  :',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','right',...
    'Position',[0.45 0.725 0.125 0.15]);

hPIVImgStart = uicontrol(...
    'Style','edit',...
    'String','1',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'Position',[0.58 0.945 0.05 0.04]);

hPIVImgEnd = uicontrol(...
    'Style','edit',...
    'String','2',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'Position',[0.58 0.895 0.05 0.04]);

hPIVImgInc = uicontrol(...
    'Style','edit',...
    'String','1',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'Position',[0.58 0.845 0.05 0.04]);

hAdjustContrastPIVRadioButton = uicontrol(...
    'Style','radiobutton',...
    'String','Adjust constrast',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'Position',[0.7 0.945 0.2 0.025]);

hUseROIPIVRadioButton = uicontrol(...
    'Style','radiobutton',...
    'String','Use ROI',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'Position',[0.7 0.895 0.2 0.025]);

hDeformROIPIVRadioButton = uicontrol(...
    'Style','radiobutton',...
    'String','Deform ROI',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'Position',[0.7 0.845 0.2 0.025]);

hStartPIV = uicontrol(...
    'Callback', @hStartPIVCallback,...
    'Style','pushbutton',...
    'String','Start',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'Position',[0.6 0.7 0.1 0.05]);

% cross correlation settings

htextPIVTitle = uicontrol(...
    'Style','text', ...
    'FontWeight', 'bold', ...
    'String','Cross-correlation settings:',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','left',...
    'Position',[0.01 0.780 0.35 0.15]);

% number passes
htextPIV3 = uicontrol(...
    'Style','text', ...
    'String','Number passes:',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','right',...
    'Position',[0.01 0.735 0.15 0.15]);
hpopupNumberPass = uicontrol(...
    'Style','popupmenu',...
    'String',{'1', '2','3','4'},...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'Position',[0.175,0.745,0.15,0.15]);

% pass 1
htextPIV1 = uicontrol(...
    'Style','text', ...
    'String','Int. Area pass 1:',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','right',...
    'Position',[0.01 0.680 0.15 0.15]);
hpopupIntAreaPass1 = uicontrol(...
    'Style','popupmenu',...
    'String',{'512', '256','128','64','32','16'},...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'Position',[0.175,0.690,0.15,0.15]);

% passe 2
htextPIV4 = uicontrol(...
    'Style','text', ...
    'String','Int. Area pass 2:',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','right',...
    'Position',[0.01 0.635 0.15 0.15]);
hpopupIntAreaPass2 = uicontrol(...
    'Style','popupmenu',...
    'String',{'512', '256','128','64','32','16'},...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'Position',[0.175,0.645,0.15,0.15]);

% passe 3
htextPIV5 = uicontrol(...
    'Style','text', ...
    'String','Int. Area pass 3:',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','right',...
    'Position',[0.01 0.590 0.15 0.15]);
hpopupIntAreaPass3 = uicontrol(...
    'Style','popupmenu',...
    'String',{'512', '256','128','64','32','16'},...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'Position',[0.175,0.600,0.15,0.15]);

% passe 4
htextPIV6 = uicontrol(...
    'Style','text', ...
    'String','Int. Area pass 4:',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','right',...
    'Position',[0.01 0.545 0.15 0.15]);
hpopupIntAreaPass4 = uicontrol(...
    'Style','popupmenu',...
    'String',{'512', '256','128','64','32','16'},...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'Position',[0.175,0.555,0.15,0.15]);

% Wind def and sub pix
htextPIV7 = uicontrol(...
    'Style','text', ...
    'String','Int. A defm:',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','right',...
    'Position',[0.01 0.490 0.15 0.15]);
hpopupWindowDeform = uicontrol(...
    'Style','popupmenu',...
    'String',{'*linear', '*spline'},...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'Position',[0.175,0.50,0.15,0.15]);

% subpixel
htextPIV8 = uicontrol(...
    'Style','text', ...
    'String','Subpixel finder:',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','right',...
    'Position',[0.01 0.445 0.15 0.15]);
hpopupSubPix = uicontrol(...
    'Style','popupmenu',...
    'String',{'3 pt Gauss', '2D Gauss'},...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'Position',[0.175,0.455,0.15,0.15]);

% filtering options
htextPIVPostTitle = uicontrol(...
    'Style','text', ...
    'FontWeight', 'bold', ...
    'String','PIV post-processing:',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','left',...
    'Position',[0.01 0.175 0.25 0.15]); %'Position',[0.01 0.385 0.25 0.15]);

% velocity filter 
hRadioVelFilter = uicontrol(...
    'Style','radiobutton',...
    'Value',1,...
    'String','Min/max filter',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','left',...
    'Position',[0.01 0.250 0.4 0.025]);
htextPIVpost1 = uicontrol(...
    'Style','text', ...
    'String','Vx min:',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','right',...
    'Position',[0.01 0.08 0.075 0.15]);
htextPIVpost2 = uicontrol(...
    'Style','text', ...
    'String','Vx max:',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','right',...
    'Position',[0.01 0.035 0.075 0.15]);
htextPIVpost3 = uicontrol(...
    'Style','text', ...
    'String','Vy min:',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','right',...
    'Position',[0.01 -0.01 0.075 0.15]);
htextPIVpost4 = uicontrol(...
    'Style','text', ...
    'String','Vy max:',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','right',...
    'Position',[0.01 -0.055 0.075 0.15]);

hPIVPostUmin = uicontrol(...
    'Style','edit',...
    'String','-10',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'Position',[0.1,0.2,0.05,0.035]); %'Position',[0.175,0.465,0.05,0.035]);
hPIVPostUmax = uicontrol(...
    'Style','edit',...
    'String','10',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'Position',[0.1,0.155,0.05,0.035]);
hPIVPostVmin = uicontrol(...
    'Style','edit',...
    'String','-10',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'Position',[0.1,0.110,0.05,0.035]);
hPIVPostVmax = uicontrol(...
    'Style','edit',...
    'String','10',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'Position',[0.1,0.065,0.05,0.035]);

% Univ Outlier
hRadioUOFilter = uicontrol(...
    'Style','radiobutton',...
    'Value',1,...
    'String','Universal outlier',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','left',...
    'Position',[0.2 0.250 0.4 0.025]);
htextPIVpost5 = uicontrol(...
    'Style','text', ...
    'String','Threshold :',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','right',...
    'Position',[0.20 0.205 0.075 0.025]);

htextPIVpost6 = uicontrol(...
    'Style','text', ...
    'String','Epsilon :',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','right',...
    'Position',[0.2 0.160 0.075 0.025]);
hPIVPostThresMed = uicontrol(...
    'Style','edit',...
    'String','2',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'Position',[0.290,0.205,0.05,0.035]);
hPIVPostEpsMed = uicontrol(...
    'Style','edit',...
    'String','0.15',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'Position',[0.290,0.160,0.05,0.035]);

% Universal Outlier
htextPIVpostKUO = uicontrol(...
    'Style','text', ...
    'String','Kernel :',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','right',...
    'Position',[0.2 0.115 0.075 0.025]);
hpopupKernelUO = uicontrol(...
    'Style','popupmenu',...
    'String',{'3x3', '5x5','7x7'},...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'Position',[0.285,0.000,0.06,0.15]);

% Std dev
hRadioStdFilter = uicontrol(...
    'Style','radiobutton',...
    'Value',1,...
    'String','Std dev filter:',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','left',...
    'Position',[0.4 0.250 0.4 0.025]);
htextPIVpost7 = uicontrol(...
    'Style','text', ...
    'String','Threshold:',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','right',...
    'Position',[0.4 0.205 0.075,0.025]);
hPIVPostThresStd = uicontrol(...
    'Style','edit',...
    'String','2',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'Position',[0.490,0.205,0.05,0.035]);

% Interpolate
hRadioInterpolate = uicontrol(...
    'Style','radiobutton',...
    'Value',1,...
    'String','Interpolate vectors:',...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'HorizontalAlignment','left',...
    'Position',[0.6 0.250 0.4 0.025]);
hpopupInterpol = uicontrol(...
    'Style','popupmenu',...
    'String',{'linear','spline','kriging','plate_0','plate_1','plate_2','plate_3','spring','Average','no interp'},...
    'Value',7,...
    'Parent', hpanelPIVSettings, ...
    'Units', 'normalized', ...
    'Position',[0.60,0.0900,0.12,0.15]);

% %% Change display setting background UI
%TecPIV_Design_DisplayFigure(hDisplayFigure, myhandles)



%% Export All Settings UI
hpanelExportAllSettings = uipanel(...
    'Title','Export settings',...
    'FontSize',8,...
    'Parent', hSecondFigure, ...
    'Units', 'normalized', ...
    'Position',[0.01 0.01 0.98 0.98]);

hExportAllPDFRadioButton = uicontrol(...
    'Style','radiobutton',...
    'String','Export as pdf',...
    'Parent', hpanelExportAllSettings, ...
    'Units', 'normalized', ...
    'Position',[0.01 0.2 0.2 0.2]);

hExportAllPNGRadioButton = uicontrol(...
    'Style','radiobutton',...
    'String','Export as png',...
    'Parent', hpanelExportAllSettings, ...
    'Units', 'normalized', ...
    'Position',[0.01 0.4 0.2 0.2]);

hExportAllAVIRadioButton = uicontrol(...
    'Style','radiobutton',...
    'String','Export as avi',...
    'Parent', hpanelExportAllSettings, ...
    'Units', 'normalized', ...
    'Position',[0.01 0.6 0.2 0.2]);

hStartExportAllButton = uicontrol(...
    'Callback', @hStartExportAll,...
    'Style','pushbutton',...
    'String','Start',...
    'Parent', hpanelExportAllSettings, ...
    'Units', 'normalized', ...
    'Position',[0.6 0.7 0.1 0.05]);


end %//folding

%





%hpanelEulerianSum.Visible='off';

%hpanelChangeBackgroundDisplaySettings.Visible='off';
%hpanelChangeVectorDisplaySettings.Visible='off';
%hpanelChangeVectorDerivativeDisplaySettings.Visible='off';
%hpanelChangeApplyDisplaySettings.Visible = 'off';

%




%% Callback functions 
% MenuItems
function hOpenMenuitemCallback(hObject,eventdata)
	warning('off','parallel:gpu:device:DeviceObjectLoaded')
      
    % Callback function run when the Open menu item is selected
    [FileName,PathName] = uigetfile('*.mat','Select the mat file');
        
    if ~isequal(FileName, 0)
        cd(PathName)
        % Delete the newly created figure
        delete(findobj(0,'type','figure'));
             
        % Load the two figures
        MainFigure;
             
        % get handle of the main figure
        hMainFigure = gcf;
            
        load(FileName, 'Data');
        fields = fieldnames(Data);
                for i=1:numel(fields)
                  myhandles.(fields{i})=Data.(fields{i});
                end   
    end
        
   guidata(hMainFigure,myhandles);
   
   if ismac
    % Code to run on Mac plaform
    fprintf('\n')
        disp('CUDA is not installed on this machine.')
    elseif isunix
    % Code to run on Linux plaform
    fprintf('\n')
    disp('CUDA is not installed on this machine.')
    elseif ispc
    % Code to run on Windows platform
    % need to reload the gpu capability
    f=fullfile(TecPivFolder ,'CUDA.mat');
    No_CUDA_msg=load(f,'No_CUDA_msg');
    [status,cmdout] = system('nvcc');

    tf = strcmp(No_CUDA_msg,cmdout)
    MaxComp=0;
    
    if tf == 1
        fprintf(' No.\n')
        fprintf('\n')
        disp('CUDA is not installed on this machine. Either there is no suitable GPU for CUDA compute, or CUDA is not installed properly.')

    else
    
        fprintf(' Yes.\n')
        fprintf('\n')
        g = gpuDevice;
    
        for ii = 1:gpuDeviceCount
            g = gpuDevice(ii);
            if g.ComputeCapability >= MaxComp
                MaxComp=g.ComputeCapability;
                Index=g.Index;
            end
        
        end
    
        gpuDevice(Index);
        fprintf('GPU Device %i has been selected\n', ...
        g.Index)
     
    end
    
    
    else
    disp('Platform not supported')
   end
        
        guidata(hMainFigure,myhandles);
end
function hNewMenuitemCallback(hNewMenuitem,eventdata)
	% Callback function run when the NewProject menu item is selected
	pathname = uigetdir('Create a new project directory - Project ID is directory name');
	cd( pathname );
	[myhandles.PathData,myhandles.ProjectID]=fileparts(pwd);
	myhandles.ProjectIDtag=[myhandles.ProjectID '.mat']; % create the mat file to save the workspace
	myhandles.NumberOfDatasets=0;
    hMainFigure.Name = ['TecPIV project: ',myhandles.ProjectID];
	guidata(hMainFigure,myhandles);
end
function hImportExpFramesMenuitemCallback(hObject,eventdata)
    % Callback function run when the Import Experiment farmes menu item is selected
    % make second figure visible
        hImportFigure.Visible = 'on';
        myhandles.ImgType='Raw';

    % make ui panel visible
        %hpanelImportFrames.Visible = 'on';
        guidata(hMainFigure,myhandles);
    end
function hImportCalibFramesMenuitemCallback(hObject,eventdata)
        % Callback function run when the Import Experiment farmes menu item is selected
        % make second figure visible
        hImportFigure.Visible = 'on';
        myhandles.ImgType='Calibration';

        % make ui panel visible
        hpanelImportFrames.Visible = 'on';
        guidata(hMainFigure,myhandles);
    end
function hExtractControlPointsMenuitemCallback(hObject,eventdata)
        hCPFigure.Visible = 'on';
        %hpanelGetControlPoints.Visible = 'on';
        guidata(hMainFigure,myhandles);
    end
function hRectifyCalibrationMenuitemCallback(hRectifyCalibrationMenuitem,eventdata)
        hSecondFigure.Visible = 'on';
        hpanelRectifySettings.Visible = 'on';
        guidata(hMainFigure,myhandles);
end
function hImportFramesCallback(hImportFrames,eventdata)
        
        warning('off', 'MATLAB:DELETE:FileNotFound') 
        
        % get camera number
        %CamFolderName=['Cam_',hCamera.String];
        %mkdir( CamFolderName);  cd( CamFolderName); % make folder Cam_1
        
        DataTypeFolderName=myhandles.ImgType; % Raw or Calibration
        mkdir( DataTypeFolderName );  % make folder
        cd('..'); % we are back in project folder
    
        ThisDataSetNumber=myhandles.NumberOfDatasets+1;
        myhandles.NumberOfDatasets=ThisDataSetNumber;
    
        ThisDataSetName=fullfile(DataTypeFolderName);
        
        myhandles.DataSets{ThisDataSetNumber,1}=ThisDataSetName;
        myhandles.DataSets{ThisDataSetNumber,2}=myhandles.PathData;
        myhandles.DataSets{ThisDataSetNumber,3}=myhandles.ProjectID;
    
        guidata(hMainFigure,myhandles);
   
        % define which format we want to import
        if hradiobuttonImportDNG.Value==1
            Iformat='.dng';
        else
            Iformat='.tif';
        end
    
        [ NumberImages,ImageWidth,ImageHeight ] = TecPIV_Import(myhandles.DataSets,ThisDataSetNumber,Iformat,myhandles.TecPivFolder);
    
        myhandles.DataSets{ThisDataSetNumber,4} = NumberImages;
        myhandles.DataSets{ThisDataSetNumber,5} = ImageWidth;
        myhandles.DataSets{ThisDataSetNumber,6} = ImageHeight;
        myhandles.DataSets{ThisDataSetNumber,7} = str2double(hTimeInc.String); % TimeIncrement
        myhandles.DataSets{ThisDataSetNumber,8} = 1; % ImageIncrement
        myhandles.DataSets{ThisDataSetNumber,9} = 1; % StartImage
        myhandles.DataSets{ThisDataSetNumber,10} = NumberImages; % EndImage
        myhandles.DataSets{ThisDataSetNumber,11} = ThisDataSetNumber; % dataset number for images associated with dataset
        myhandles.DataSets{ThisDataSetNumber,12} = 1; % image scale %ThisDataSetNumber;
        myhandles.DataSets{ThisDataSetNumber,13} = 'pix'; % whether image is calibrated or not
        myhandles.DataSets{ThisDataSetNumber,14} = 'pix'; % physical para for scale
        myhandles.DataSets{ThisDataSetNumber,15} = 'image'; % type of dataset
        myhandles.DataSets{ThisDataSetNumber,16} = 0;
        myhandles.DataSets{ThisDataSetNumber,17} = [];
    
        guidata(hMainFigure,myhandles);
    
        % Display the first image
        TecPIV_Display_2(hPlotAxes,1,myhandles.DataSets,ThisDataSetNumber,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
        
        %Frame=1;
        %FramePath=fullfile(myhandles.PathData,myhandles.ProjectID,ThisDataSetName,['IMG_' num2str(Frame) '.tif']);
        %I0=imread(FramePath);
        %TecPIV_Display(myhandles.TecPivFolder,I0,hPlotAxes,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
    
        % Update the player
        [ ImgNumber, Time, Min, Max, SliderStep , Value] = TecPIV_UpdatePlayer('FrameNum','1',myhandles.DataSets,ThisDataSetNumber); 
        myhandles.ImgNumber=ImgNumber;
        myhandles.Time=Time;
    
        if NumberImages == 1
            hslider.Visible='off';
            hImgNumber.Visible='off';
            hImgTime.Visible='off';
        else
            
            hslider.Visible='on';
            hImgNumber.String=num2str(myhandles.ImgNumber);
            hImgTime.String=num2str(myhandles.Time);
            
            SliderStep=1/(ceil(NumberImages/1)-1);
            SliderStep=[SliderStep SliderStep];
            SliderMax=ceil(NumberImages/1);
            SliderMin=1;
            
            
            hslider.SliderStep=SliderStep;
            hslider.Max=Max;
            hslider.Min=Min;
            hslider.Value=Value;

        end
    
        % Add entry to the source selector
        NewEntryName=myhandles.DataSets{ThisDataSetNumber,1};
        myhandles.entries = hpopupSourceSelector.String;
        myhandles.entries = [myhandles.entries; NewEntryName];
        hpopupSourceSelector.String = myhandles.entries;
        myhandles.SelectedEntry=ThisDataSetNumber;
        hpopupSourceSelector.Value=ThisDataSetNumber; % Select new Entry
    
        hpanelImportFrames.Visible = 'off';
        hSecondFigure.Visible = 'off';
        set(hMainFigure, 'pointer','arrow')
    
        guidata(hMainFigure,myhandles);
    
end
function hControlPointsExtractionStartCallback(hStartControlPointsExtraction,eventdata)
        % get the parameters for the calibration board from GUI
        myhandles.n_sq_x=str2double(hNbSqX.String);
        myhandles.n_sq_y=str2double(hNbSqY.String);
        myhandles.dX=str2double(hSizeSqX.String);
        myhandles.dY=str2double(hSizeSqY.String);
        myhandles.PhysU=hPhysU.String;

        % Find the dataset number for the dataset Calibration 
        String=fullfile('Calibration');
        index = find(strcmp(myhandles.DataSets, String));
        DataSetNumber=index(1);

        myhandles.DataSets;
        guidata(hMainFigure,myhandles);
        
        TecPIV_ExtractGCP(myhandles.TecPivFolder,myhandles.DataSets,DataSetNumber,hPlotAxes,myhandles.RawCpt,myhandles.n_sq_x,myhandles.n_sq_y,myhandles.dX,myhandles.dY,myhandles.VectorField,myhandles.Derivative );
        hpanelGetControlPoints.Visible = 'off';
        hSecondFigure.Visible = 'off';
        
        set(hMainFigure, 'pointer','arrow')
        guidata(hMainFigure,myhandles);
    
end
function hStartRectifyCalibCallback(hStartRectifyCalib,eventdata)
        %get the parameters
        FrameNum=str2double(hRectFrameNum.String);
        Order=str2double(hOrderPoly.String);
        RectMethod=hpopupRectMethodSelector.Value;
        
        %find the dataset
        String=fullfile('Calibration');
        index = find(strcmp(myhandles.DataSets, String));
        DataSetNumber=index(1);
        
        [ImScale,RectFn,tx,wintx,winty] = TecPIV_Rectify(myhandles.DataSets,DataSetNumber,FrameNum,RectMethod,Order);
        myhandles.ImScale=ImScale;
        myhandles.RectFn=RectFn;     

        DataFolder=myhandles.DataSets{DataSetNumber,1};
        NewDataFolder=fullfile(DataFolder,'Rectified');

        % Add entry to the source selector for the calibrated calibration
        % image
        ThisDataSetNumber=myhandles.NumberOfDatasets+1;
        myhandles.NumberOfDatasets=ThisDataSetNumber;

        myhandles.DataSets{ThisDataSetNumber,1}=NewDataFolder;
        myhandles.DataSets{ThisDataSetNumber,2}=myhandles.PathData;
        myhandles.DataSets{ThisDataSetNumber,3}=myhandles.ProjectID;
        
        myhandles.DataSets{ThisDataSetNumber,11}=ThisDataSetNumber;
        myhandles.DataSets{ThisDataSetNumber,15}='image';

        guidata(hMainFigure);
    
        % Display the first image
        TecPIV_Display_2(hPlotAxes,1,myhandles.DataSets,ThisDataSetNumber,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
        
%         Frame=1;
%         FramePath=fullfile(myhandles.PathData,myhandles.ProjectID,NewDataFolder,['IMG_' num2str(Frame) '.tif']);
%         I0=imread(FramePath);
%         TecPIV_Display(myhandles.TecPivFolder,I0,hPlotAxes,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
        
        
        n_sq_x=str2double(hNbSqX.String);
        n_sq_y=str2double(hNbSqY.String);
        dX=str2double(hSizeSqX.String);
        dY=str2double(hSizeSqY.String);
        
           if RectMethod == 3
               STEP=1/2;
           else
               STEP=1/1;
           end
           
        Frame=1;
        FramePath=fullfile(myhandles.PathData,myhandles.ProjectID,NewDataFolder,['IMG_' num2str(Frame) '.tif']);
        I0=imread(FramePath);
        [SDX, SDY,CropImage,rect] = TecPIV_CheckCalib(I0,...
            hPlotAxes,...
            myhandles.RawCpt,...
            tx,...
            wintx,...
            winty,...
            n_sq_x,...
            n_sq_y,...
            dX,...
            dY,...
            myhandles.PathData,...
            myhandles.ProjectID,...
            NewDataFolder,...
            ImScale,...
            myhandles.VectorField,...
            myhandles.Derivative,...
            STEP);
        
        myhandles.CropImage=CropImage;
        myhandles.CropRect=rect;
        
        if RectMethod == 3
            STEP=2/2;
            
            % now we do second step => polynomial
            [RectFn2] = TecPIV_Rectify_second_step(myhandles.DataSets,ThisDataSetNumber,1,2,Order,ImScale);
            myhandles.RectFn2=RectFn2;

            % Display the first image
            TecPIV_Display_2(hPlotAxes,1,myhandles.DataSets,ThisDataSetNumber,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
            
%             Frame=1;
%             NewDataFolder=fullfile(NewDataFolder,'Rectified2');
%             FramePath=fullfile(myhandles.PathData,myhandles.ProjectID,NewDataFolder,['IMG_' num2str(Frame) '.tif']);
%             I0=imread(FramePath);
%             TecPIV_Display(myhandles.TecPivFolder,I0,hPlotAxes,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
            Frame=1;
            FramePath=fullfile(myhandles.PathData,myhandles.ProjectID,NewDataFolder,['IMG_' num2str(Frame) '.tif']);
            I0=imread(FramePath);
            [SDX, SDY,CropImage,rect] = TecPIV_CheckCalib(I0,...
                hPlotAxes,...
                myhandles.RawCpt,...
                tx,...
                wintx,...
                winty,...
                n_sq_x,...
                n_sq_y,...
                dX,...
                dY,...
                myhandles.PathData,...
                myhandles.ProjectID,...
                NewDataFolder,...
                ImScale,...
                myhandles.VectorField,...
                myhandles.Derivative,...
                STEP);
            myhandles.CropImage=CropImage;
            myhandles.CropRect=rect;
            
            
        else
        end
        
        SizeI0=size(I0);
        ImageHeight=SizeI0(1,1);
        ImageWidth=SizeI0(1,2);

        myhandles.DataSets{ThisDataSetNumber,4} = 1;
        myhandles.DataSets{ThisDataSetNumber,5} = ImageWidth;%%
        myhandles.DataSets{ThisDataSetNumber,6} = ImageHeight;%%%
        myhandles.DataSets{ThisDataSetNumber,7} = 1; % dt (time inc)
        myhandles.DataSets{ThisDataSetNumber,8} = 1; % image increment
        myhandles.DataSets{ThisDataSetNumber,9} = 1; % start number
        myhandles.DataSets{ThisDataSetNumber,10} = 1; % end number       
        myhandles.DataSets{ThisDataSetNumber,11} = ThisDataSetNumber; % number of the dataset containing the images (self if this dataset is image)
        myhandles.DataSets{ThisDataSetNumber,12} = myhandles.ImScale;
        myhandles.DataSets{ThisDataSetNumber,13} = 'phys'; % toggle scalebar on/off
        myhandles.DataSets{ThisDataSetNumber,14} = myhandles.PhysU;
        myhandles.DataSets{ThisDataSetNumber,15} = 'image';
        myhandles.DataSets{ThisDataSetNumber,16} = 0; % no mask
        myhandles.DataSets{ThisDataSetNumber,17} = []; % empty cell because no mask
        
        
        myhandles.RawCpt{1,4}=myhandles.DataSets{ThisDataSetNumber,13}; %togle scalebar on/off
        myhandles.RawCpt{1,5}=myhandles.DataSets{ThisDataSetNumber,12}; % ImScale
        myhandles.RawCpt{1,6}= myhandles.DataSets{ThisDataSetNumber,14}; % physical unit
        
        NewEntryName=myhandles.DataSets{ThisDataSetNumber,1};
        entries= hpopupSourceSelector.String;
        entries = [entries; NewEntryName];
        hpopupSourceSelector.String = entries;
        hpopupSourceSelector.Value=ThisDataSetNumber; % Select new Entry
        
        hpanelRectifySettings.Visible = 'off';
        hSecondFigure.Visible = 'off';
        
        set(hMainFigure, 'pointer','arrow')
        guidata(hMainFigure,myhandles);
end
function hDisplaySettingsMenuitemCallback(hDisplaySettingsMenuitem,eventdata)
        % make second figure visible
        hDisplayFigure.Visible = 'on';

        % make ui panel visible
        hpanelChangeBackgroundDisplaySettings.Visible = 'on';
        hpanelChangeVectorDisplaySettings.Visible = 'on';
        hpanelChangeVectorDerivativeDisplaySettings.Visible = 'on';
        %hpanelChangeApplyDisplaySettings.Visible = 'on';
        % place the variables in GUI panel
        % Background image
        % list of palette names
        ListPopUp = hpopupBackgroundCPTSelector.String;
        LenList=length(ListPopUp);
        Choice= cellstr(myhandles.RawCpt{1,1}); % selected palette name
        
        for i=1:LenList
            Palette=char(ListPopUp(i));
            test = strcmp(Palette,Choice);
            if test == 1
                 hpopupBackgroundCPTSelector.Value=i; % actually selected palette
            end
        end
        
        hMinBackgroundColorPalette.String = myhandles.RawCpt{1,2};
        hMaxBackgroundColorPalette.String = myhandles.RawCpt{1,3};
        
        % vectors
        hDisplayVectorRadioButton.Value = myhandles.VectorField{1,1};
        
        hpopupVecColorSelector.Value=8; % black default
        
        temp = myhandles.VectorField{1,2};
        hpopupVecGridFactorSelector.Value = temp(1,1);
        
        DisplayMode=myhandles.VectorField{1,3};
        
        switch DisplayMode
            case 'max'
                hpopupVecColorScalingModeSelector.Value = 2;
            case 'mean'
                hpopupVecColorScalingModeSelector.Value = 1;
            case 'manual'
                hpopupVecColorScalingModeSelector.Value = 3;
        end
        
        % derivatives
        hDisplayVectorDerivativeRadioButton.Value = myhandles.Derivative{1,1}; % display yes/no
        
        ListPopUpDeriv =  hpopupVecDerivativeTypeSelector.String;
        LenList=length(ListPopUpDeriv);
        Choice= myhandles.Derivative{1,2}; % selected programmatically
        hpopupVecDerivativeTypeSelector.Value=1; % default
        
        for i=1:LenList
           
            DerivType=char(ListPopUpDeriv(i));
            test = strcmp(DerivType,Choice);
            if test == 1
                 hpopupVecDerivativeTypeSelector.Value=i; % actually selected derivative
            end
        end
        
        hpopupVecDerivativeDisplayRangeSelector.Value=myhandles.Derivative{1,3}; %1 = minmax; 2 = +- max; 3 = custom
       
        ListCPTDeriv=hpopupVecDerivativeCPTSelector.String; % list of CPT for derivative
        LenList=length(ListCPTDeriv);
        Choice= myhandles.Derivative{1,4}; % selected derivative cpt
        for i=1:LenList
            DerivCPT=strtrim(char(ListCPTDeriv(i)));
            test = strcmp(DerivCPT,Choice);
            if test == 1
                hpopupVecDerivativeCPTSelector.Value = i;
            end
        end
        
        hMinVecDerivativeColorPalette.String = num2str(myhandles.Derivative{1,5});
        hMaxVecDerivativeColorPalette.String = num2str(myhandles.Derivative{1,6});
        
        guidata(hMainFigure,myhandles);
    
end
function hExportThisFrameMenuitemCallback (hExportThisFrameMenuitem,eventdata)
    
        SourceNum=hpopupSourceSelector.Value;
        SourceName=char(hpopupSourceSelector.String(SourceNum))
        
        ThisDataSetNumber=hpopupSourceSelector.Value; %get  the selected dataset
        CurrentFrame=str2num(hImgNumber.String); %get the selected frame
        DatasetFolder = myhandles.DataSets{ThisDataSetNumber,1};
        PathData = myhandles.DataSets{ThisDataSetNumber,2};
        ProjectID = myhandles.DataSets{ThisDataSetNumber,3};
        NumberImages = myhandles.DataSets{ThisDataSetNumber,4};
        ImageWidth = myhandles.DataSets{ThisDataSetNumber,5};%%
        ImageHeight = myhandles.DataSets{ThisDataSetNumber,6};%%%
        TimeInc = myhandles.DataSets{ThisDataSetNumber,7};
        ImageInc = myhandles.DataSets{ThisDataSetNumber,8};
        
        % check if export folder exist
        NewExportFolderName = 0;
        k=1;
        while NewExportFolderName == 0
        ExportFolderName = fullfile(PathData,ProjectID,['Export_',num2str(k)]);
            if ~exist(ExportFolderName,'dir')
                NewExportFolderName = 1;
                mkdir(ExportFolderName);
            else
                k=k+1;
            end
        end
        
        % create new figure
        h = figure();
        newaxes = axes('DataAspectRatio',[1 1 1]);
        % Display image
        TecPIV_Display_2(newaxes,CurrentFrame,myhandles.DataSets,ThisDataSetNumber,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
        
        Colormap = myhandles.Derivative{1,4}
        ListMATLABCPT={'parula','jet','hsv','hot','cool','spring','summer',...
                'autumn','winter','gray','bone','copper','pink','lines',...
                'colorcube','prism','flag','white'};
        Lia = ismember(Colormap,ListMATLABCPT); % check if colormap is in the list of MATLAB colormaps
                
        if Lia == 1 % it is a MATLAB colormap
                    RGB=Colormap;
                    colormap(newaxes,RGB);
        else % it is a cutom colormap loaded from the folder
                    load(fullfile(myhandles.TecPivFolder,'toolbox','colormaps',Colormap));
                    colormap(newaxes,RGB);
        end 
        
        
        % Set up the paper size / position
        set(h, 'PaperUnits', 'centimeters');
        set(h, 'PaperSize', [21 29.7]);    % Set to final desired size here as well as 2 lines below
        set(h, 'PaperPositionMode', 'manual');
        set(h, 'PaperPosition', [0 0 21 29.7]);
        
        switch SourceName       
            case 'Raw'
                filename=strcat('Raw-IMG',hImgNumber.String,'.pdf')  
            case 'Calibration'
                filename=strcat('Calibration-IMG',hImgNumber.String,'.pdf')  
            case (fullfile('Calibration','Rectified'))
                filename=strcat('Calibration-Rectified-IMG',hImgNumber.String,'.pdf')
            case (fullfile('Raw','Rectified'))
                filename=strcat('Raw-Rectified-IMG',hImgNumber.String,'.pdf')
            case (fullfile('Raw','Rectified','Vectors'))
                filename=strcat('Raw-Rectified-Vectors-IMG_',hImgNumber.String,'.pdf')
            case (fullfile('Raw','Vectors'))
                filename=strcat('Raw-Vectors-IMG_',hImgNumber.String,'.pdf')
        end

       print(h,'-painters','-dpdf',fullfile(ExportFolderName,filename))  
       close(h)
       guidata(hMainFigure,myhandles);
        
        
%         c1 = fullfile('Raw');
%         c2 = fullfile('Calibration');
%         c3 = fullfile('Calibration','Rectified');
%         c4 = fullfile('Raw','Rectified');
%         c5 = fullfile('Raw','Rectified','Vectors');
%         c6 = fullfile('Raw','Vectors');
%         
%         choices = {c1,c2,c3,c4,c5,c6};
%         test = strcmp(SourceName,choices);
%         
%         if test == 1
%             String = c1;
%             index = find(strcmp(myhandles.DataSets, String));
%             BackgroundImageDataSetNumber=index(1); 
%             BackgroundImageDatasetFolder = myhandles.DataSets{BackgroundImageDataSetNumber,1};
%         
%             FramePath=fullfile(PathData,ProjectID,BackgroundImageDatasetFolder,['IMG_' num2str(CurrentFrame) '.tif']);
%             I0=imread(FramePath);
%             
%             filename=strcat('Raw-IMG',hImgNumber.String,'.pdf')
%             
%         elseif test == 2
%             String = c2;
%             index = find(strcmp(myhandles.DataSets, String));
%             BackgroundImageDataSetNumber=index(1); 
%             BackgroundImageDatasetFolder = myhandles.DataSets{BackgroundImageDataSetNumber,1};
%         
%             FramePath=fullfile(PathData,ProjectID,BackgroundImageDatasetFolder,['IMG_' num2str(CurrentFrame) '.tif']);
%             I0=imread(FramePath);
%             
%             filename=strcat('Calibration-IMG',hImgNumber.String,'.pdf')
%             
%         elseif test == 3
%             String = c3;
%             index = find(strcmp(myhandles.DataSets, String));
%             BackgroundImageDataSetNumber=index(1); 
%             BackgroundImageDatasetFolder = myhandles.DataSets{BackgroundImageDataSetNumber,1};
%         
%             FramePath=fullfile(PathData,ProjectID,BackgroundImageDatasetFolder,['IMG_' num2str(CurrentFrame) '.tif']);
%             I0=imread(FramePath);
%             
%             filename=strcat('Calibration-Rectified-IMG',hImgNumber.String,'.pdf')
%             
%         elseif test == 4
%             String = c4;
%             index = find(strcmp(myhandles.DataSets, String));
%             BackgroundImageDataSetNumber=index(1); 
%             BackgroundImageDatasetFolder = myhandles.DataSets{BackgroundImageDataSetNumber,1};
%         
%             FramePath=fullfile(PathData,ProjectID,BackgroundImageDatasetFolder,['IMG_' num2str(CurrentFrame) '.tif']);
%             I0=imread(FramePath);
%             
%             filename=strcat('Raw-Rectified-IMG',hImgNumber.String,'.pdf')
%             
%         elseif test == 5
%             String = c4; % image is raw rectified
%             index = find(strcmp(myhandles.DataSets, String));
%             BackgroundImageDataSetNumber=index(1); 
%             BackgroundImageDatasetFolder = myhandles.DataSets{BackgroundImageDataSetNumber,1};
%         
%             FramePath=fullfile(PathData,ProjectID,BackgroundImageDatasetFolder,['IMG_' num2str(CurrentFrame) '.tif']);
%             I0=imread(FramePath);
%             
%             % get vector field
%             load(fullfile(PathData,ProjectID,BackgroundImageDatasetFolder,'Vectors', ['Vector_' num2str(CurrentFrame) '.mat']),'X','Y','U','V');  
%             myhandles.VectorField{1,5} = X; %cell2mat(X); 
%             myhandles.VectorField{1,6} = Y; %cell2mat(Y); 
%             myhandles.VectorField{1,7} = U; %cell2mat(U); 
%             myhandles.VectorField{1,8} = V; %cell2mat(V);
%             
%             filename=strcat('Raw-Rectified-Vectors-IMG',hImgNumber.String,'.pdf')
%             
%         else
%             String = c1; % image is raw
%             index = find(strcmp(myhandles.DataSets, String));
%             BackgroundImageDataSetNumber=index(1); 
%             BackgroundImageDatasetFolder = myhandles.DataSets{BackgroundImageDataSetNumber,1}
%         
%             FramePath=fullfile(PathData,ProjectID,BackgroundImageDatasetFolder,['IMG_' num2str(CurrentFrame) '.tif']);
%             I0=imread(FramePath);
%             
%             % get vector field
%             load(fullfile(PathData,ProjectID,BackgroundImageDatasetFolder,'Vectors', ['Vector_' num2str(CurrentFrame) '.mat']),'X','Y','U','V');  
%             myhandles.VectorField{1,5} = X; %cell2mat(X); 
%             myhandles.VectorField{1,6} = Y; %cell2mat(Y); 
%             myhandles.VectorField{1,7} = U; %cell2mat(U); 
%             myhandles.VectorField{1,8} = V; %cell2mat(V); 
%             
%             filename=strcat('Raw-Vectors-IMG',hImgNumber.String,'.pdf')
%             
%         end
%             
%         
%         
% %         switch SourceName
% %             case 'Raw'
% %                 
% %                 String=fullfile('Raw');
% %                 index = find(strcmp(myhandles.DataSets, String));
% %                 BackgroundImageDataSetNumber=index(1); 
% %                 BackgroundImageDatasetFolder = myhandles.DataSets{BackgroundImageDataSetNumber,1};
% %         
% %                 FramePath=fullfile(PathData,ProjectID,BackgroundImageDatasetFolder,['IMG_' num2str(CurrentFrame) '.tif']);
% %                 I0=imread(FramePath);
% %                 
% %             case 'Calibration'
% %                 
% %                 String=fullfile('Calibration');
% %                 index = find(strcmp(myhandles.DataSets, String));
% %                 BackgroundImageDataSetNumber=index(1); 
% %                 BackgroundImageDatasetFolder = myhandles.DataSets{BackgroundImageDataSetNumber,1};
% %         
% %                 FramePath=fullfile(PathData,ProjectID,BackgroundImageDatasetFolder,['IMG_' num2str(CurrentFrame) '.tif']);
% %                 I0=imread(FramePath);
% %                 
% %            case 'Calibration\Rectified'
% %                 
% %                 String=fullfile('Calibration','Rectified');
% %                 index = find(strcmp(myhandles.DataSets, String));
% %                 BackgroundImageDataSetNumber=index(1); 
% %                 BackgroundImageDatasetFolder = myhandles.DataSets{BackgroundImageDataSetNumber,1};
% %         
% %                 FramePath=fullfile(PathData,ProjectID,BackgroundImageDatasetFolder,['IMG_' num2str(CurrentFrame) '.tif']);
% %                 I0=imread(FramePath);
% %                 
% %             case 'Raw\Rectified'
% %                 
% %                 String=fullfile('Raw','Rectified');
% %                 index = find(strcmp(myhandles.DataSets, String));
% %                 BackgroundImageDataSetNumber=index(1); 
% %                 BackgroundImageDatasetFolder = myhandles.DataSets{BackgroundImageDataSetNumber,1};
% %         
% %                 FramePath=fullfile(PathData,ProjectID,BackgroundImageDatasetFolder,['IMG_' num2str(CurrentFrame) '.tif']);
% %                 I0=imread(FramePath);
% %                 
% %             case 'Raw\Rectified\Vectors'
% %                 
% %                 String=fullfile('Raw','Rectified');
% %                 index = find(strcmp(myhandles.DataSets, String));
% %                 BackgroundImageDataSetNumber=index(1); 
% %                 BackgroundImageDatasetFolder = myhandles.DataSets{BackgroundImageDataSetNumber,1};
% %         
% %                 FramePath=fullfile(PathData,ProjectID,BackgroundImageDatasetFolder,['IMG_' num2str(CurrentFrame) '.tif']);
% %                 I0=imread(FramePath);
% %                 
% %                 % get vector field
% %                 load(['Raw\Rectified\Vectors\Vector_' num2str(CurrentFrame) '.mat'],'X','Y','U','V');  
% %                                
% %                 myhandles.VectorField{1,5} = X; %cell2mat(X); 
% %                 myhandles.VectorField{1,6} = Y; %cell2mat(Y); 
% %                 myhandles.VectorField{1,7} = U; %cell2mat(U); 
% %                 myhandles.VectorField{1,8} = V; %cell2mat(V); 
% %             
% %         end
% 
%                 % create new figure
%                 h = figure();
%                 newaxes = axes('DataAspectRatio',[1 1 1]);
%                 % Display image
%                 TecPIV_Display_2(newaxes,CurrentFrame,myhandles.DataSets,ThisDataSetNumber,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
%                 
%                 %TecPIV_Display(myhandles.TecPivFolder,I0,newaxes,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
%                 set(h,'Visible', 'off'); 
%         
%                 % Set up the paper size / position
%                 set(h, 'PaperUnits', 'centimeters');
%                 set(h, 'PaperSize', [21 29.7]);    % Set to final desired size here as well as 2 lines below
%                 set(h, 'PaperPositionMode', 'manual');
%                 set(h, 'PaperPosition', [0 0 21 29.7]);
%         
% %         switch SourceName
% %             
% %             case 'Raw'
% %                 filename=strcat('Raw-IMG',hImgNumber.String,'.pdf')  
% %             case 'Calibration'
% %                 filename=strcat('Calibration-IMG',hImgNumber.String,'.pdf')  
% %             case 'Calibration\Rectified'
% %                 filename=strcat('Calibration-Rectified-IMG',hImgNumber.String,'.pdf')
% %             case 'Raw\Rectified'
% %                 filename=strcat('Raw-Rectified-IMG',hImgNumber.String,'.pdf')
% %             case 'Raw\Rectified\Vectors'
% %                 filename=strcat('Raw-Rectified-Vectors-IMG',hImgNumber.String,'.pdf')
% % 
% %         end
% 
%        
%         
%        print(h,'-painters','-dpdf',filename)
%  
%         close(h)
%         guidata(hMainFigure,myhandles);
end
function hUndeformMenuitemCallback(hUndeformMenuitem,eventdata)
        % make second figure visible
        hSecondFigure.Visible = 'on';
        hpanelUndeform.Visible='on'; 
        guidata(hMainFigure,myhandles);
end
function hStartUndeformCallback (hStartUndeform,eventdata)
        
        DoRotate=hradiobuttonRotateImages.Value;
        Rotation=str2num(htextDefineAngleRotation.String);
        ImageInc=1; % process all raw images
        StartNumber = 1; % start with the first
       
        
        tic;
        %find the dataset
        String=fullfile('Raw');
        index = find(strcmp(myhandles.DataSets, String));
        DataSetNumber=index(1); 

        NumberImages=myhandles.DataSets{DataSetNumber,4};
        EndNumber = NumberImages; % since Start=1 and Inc = 1
        
        DataFolder=myhandles.DataSets{DataSetNumber,1}; 

        cd( DataFolder ); 
        mkdir( 'Rectified' );  
        cd(fullfile(myhandles.PathData,myhandles.ProjectID));
        NewDataFolder=fullfile(DataFolder,'Rectified');
               
        PathData=myhandles.PathData;
        PROJECTID=myhandles.ProjectID;
        RECTFN=myhandles.RectFn;
        
        RectMethod=hpopupRectMethodSelector.Value;
        if RectMethod == 3
                RECTFN2=myhandles.RectFn2;
        else
            RECTFN2=0;
        end
        
        CheckVariable = isfield(myhandles ,'CropImage');
        
        if CheckVariable == 1 % CropImage do exist in myhandles
            %disp('variable crop image exist')
            CROPIMAGE=myhandles.CropImage;
            CROPRECT=myhandles.CropRect;
        else
            %disp('variable crop image does not exist')
            CROPIMAGE=0;
            CROPRECT=0;
        end
        
       obj = ProgressBar(NumberImages, ...
           'IsParallel', true, ...
           'WorkerDirectory', pwd, ...
           'Title', 'Rectifying images' ...
           );

        % ALWAYS CALL THE SETUP() METHOD FIRST!!!
        obj.setup([], [], []);
        
        parfor i=1:NumberImages 
            FramePath=fullfile(PathData,PROJECTID,DataFolder,['IMG_' num2str(i) '.tif']);
            I=imread(FramePath);
            
            if RectMethod == 1 || RectMethod == 2
                I1 = imtransform(I,RECTFN,'bilinear','FillValues', Inf);
                if CROPIMAGE == 1
                    I1=imcrop(I1,CROPRECT);
                end
                if DoRotate == 1
                    I1 = imrotate(I1, Rotation, 'loose', 'bilinear');
                end
                
            else
                I1 = imtransform(I,RECTFN,'bilinear','FillValues', Inf);
                I2 = imtransform(I1,RECTFN2,'bilinear','FillValues', Inf);
                if CROPIMAGE == 1
                    I1=imcrop(I2,CROPRECT); 
                end
                if DoRotate == 1
                    I1 = imrotate(I1, Rotation, 'loose', 'bilinear');
                end
            end
            
            % save rectified image as 16-bit tiff
            cd(NewDataFolder);
            name=strcat('IMG_',num2str(i),'.tif');
            imwrite(I1,name,'tiff');

            %go back to project folder
            cd(fullfile(PathData,PROJECTID));
            
        updateParallel([], pwd);
        end
        obj.release();
        
        clear I;
        clear I0;
        UndistordTime=toc;
        message=sprintf('%0.0f images processed in %0.2f s',NumberImages,UndistordTime); 
        disp(message)
        
        hpanelUndeform.Visible='off';
        hSecondFigure.Visible = 'off';
        set(hMainFigure, 'pointer','arrow')
        
%         %% Display the first image
%         Frame=1;
        FramePath=fullfile(PathData,PROJECTID,NewDataFolder,['IMG_' num2str(StartNumber) '.tif']);
        I0=imread(FramePath);
%         TecPIV_Display(myhandles.TecPivFolder,I0,hPlotAxes,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
%         
        SizeI0=size(I0);
        ImageHeight=SizeI0(1,1);
        ImageWidth=SizeI0(1,2);
        
        %% Add entry to the source selector
        ThisDataSetNumber=myhandles.NumberOfDatasets+1; % number for the dataset with rectified images of the model
        myhandles.NumberOfDatasets=ThisDataSetNumber;

        myhandles.DataSets{ThisDataSetNumber,1} = NewDataFolder;
        myhandles.DataSets{ThisDataSetNumber,2} = myhandles.PathData;
        myhandles.DataSets{ThisDataSetNumber,3} = myhandles.ProjectID;
        myhandles.DataSets{ThisDataSetNumber,4} = NumberImages;
        myhandles.DataSets{ThisDataSetNumber,5} = ImageWidth;%%
        myhandles.DataSets{ThisDataSetNumber,6} = ImageHeight;%%%
        myhandles.DataSets{ThisDataSetNumber,7} = myhandles.DataSets{DataSetNumber,7}; % same time inc as originial dataset (for now)
        myhandles.DataSets{ThisDataSetNumber,8} = ImageInc;
        myhandles.DataSets{ThisDataSetNumber,9} = StartNumber;
        myhandles.DataSets{ThisDataSetNumber,10} = EndNumber;
        myhandles.DataSets{ThisDataSetNumber,11} = ThisDataSetNumber; % number of dataset with images associated . For images. Image source is the same dataset (not true for vectors and cumulative vectors).
        myhandles.DataSets{ThisDataSetNumber,12} = myhandles.ImScale;
        myhandles.DataSets{ThisDataSetNumber,13} = 'phys';
        myhandles.DataSets{ThisDataSetNumber,14} = myhandles.PhysU;
        myhandles.DataSets{ThisDataSetNumber,15} = 'image';
        myhandles.DataSets{ThisDataSetNumber,16} = 0; % no mask
        myhandles.DataSets{ThisDataSetNumber,17} = []; % empty cell because no mask
        
        
        % Display image
        TecPIV_Display_2(hPlotAxes,StartNumber,myhandles.DataSets,ThisDataSetNumber,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
        
        
        TimeInc = myhandles.DataSets{ThisDataSetNumber,7};
        NewEntryName=myhandles.DataSets{ThisDataSetNumber,1};
        entries= hpopupSourceSelector.String;
        entries = [entries; NewEntryName];
        hpopupSourceSelector.String = entries;
        hpopupSourceSelector.Value=ThisDataSetNumber; % Select new Entry
                
        
        %% Update the player
        Frame = StartNumber;
        
        if NumberImages == 1 % no need if only one image
            hslider.Visible='off';
            hImgNumber.Visible='off';
            hImgTime.Visible='off';
        else
            hslider.Visible='on';
            hImgNumber.Visible='on';
            hImgTime.Visible='on';
            
            SliderStep=1/(ceil(NumberImages/ImageInc)-1);
            SliderStep=[SliderStep SliderStep];
            SliderMax=ceil(NumberImages/ImageInc);
            SliderMin=1;
            
            hslider.SliderStep=SliderStep;
            hslider.Max=SliderMax;
            hslider.Min=SliderMin;
            hslider.Value=Frame;
            
            hImgNumber.String=num2str(Frame);
            CurrentTime = (Frame-StartNumber)*TimeInc; % calculate the time
            hImgTime.String=num2str(CurrentTime); % and time
            
        end
        
        %%
        guidata(hMainFigure,myhandles);
end
function hContrastMenuCallback (hContrastMenu, eventdata)
        % make second figure visible
        hSecondFigure.Visible = 'on';
        hpanelContrast.Visible='on'; 
        guidata(hMainFigure,myhandles);
end
function hApplyContrastCallback (hApplyContrast,eventdata)
        
        ThisDataSetNumber=hpopupSourceSelector.Value; %get  the selected dataset
        CurrentFrame=str2num(hImgNumber.String); %get the selected frame
        DatasetFolder = myhandles.DataSets{ThisDataSetNumber,1};
        PathData = myhandles.DataSets{ThisDataSetNumber,2};
        ProjectID = myhandles.DataSets{ThisDataSetNumber,3};
        NumberImages = myhandles.DataSets{ThisDataSetNumber,4};
        ImageWidth = myhandles.DataSets{ThisDataSetNumber,5};%%
        ImageHeight = myhandles.DataSets{ThisDataSetNumber,6};%%%
        TimeInc = myhandles.DataSets{ThisDataSetNumber,7};
                
        FramePath=fullfile(PathData,ProjectID,DatasetFolder,['IMG_' hImgNumber.String '.tif']);
        I0=imread(FramePath);
% not sure that the function below is really useful. I don't remember using it at all.
% let see if we can remove.

%         if hUseMaskRadioButton.Value == 1 % use mask
%             CheckExist=exist('ROI','var');
%             
%             if CheckExist == 1 % if mask exist
%                 % ask user if re-use existing mask
%                 % Construct a questdlg with 2 options
%                 ReUseMask = questdlg('A mask already exist. Would you like to re-use it?', ...
%                 'Yes','No, create a new one');
%                 % Handle response
%                 switch ReUseMask
%                     case 'Yes'
%                     case 'No, create a new one'
%                         uiwait(msgbox('Using the mouse, specify the region by selecting vertices of the polygon. You can move or resize the polygon using the mouse. Create the mask by double-clicking, or by right-clicking inside the region and selecting Create mask from the context menu.'))
%                         Mask = roipoly(I0);
%                         ROI = uint16(Mask);
%                 end 
%             else % ROI does not exist
%                 uiwait(msgbox('Using the mouse, specify the region by selecting vertices of the polygon. You can move or resize the polygon using the mouse. Create the mask by double-clicking, or by right-clicking inside the region and selecting Create mask from the context menu.'))
%                 Mask = roipoly(I0);
%                 ROI = uint16(Mask); 
%                 myhandles.ROI=ROI;
%             end
%    
%         I0=I0.*ROI;
%         else
%         I0=I0; 
%         end
        
        if hInverseImageRadioButton.Value == 1
            message=sprintf('Inverting Image ...');
            disp(message)
            I1=imcomplement(I0);
        else
            I1=I0;
        end
        
        if hSubtractBackgroundRadioButton.Value == 1
            message=sprintf('Subtracting Background ...');
            disp(message)
            GaussianWindowSize = str2num(hGaussianWindowSizeText.String);
            GaussianSigma = str2num(hGaussianSigmaText.String);
            hfilter=fspecial('gaussian', GaussianWindowSize, GaussianSigma);
            background=imfilter(I1,hfilter,'replicate');
            I1 = I1 - background;
            I1 = imadjust(I1);
        end
        
        if hNormalizeIntensityRadioButton.Value == 1
            message=sprintf('Normalizing intensity ...');
            disp(message)
            I1 = im2double(I1); % convert to double
            I1 = (I1 - mean(I1(:))) / std(I1(:)); %create normalized
            I1 = im2uint16(I1); %// Convert to uint8

        end
        
        % create a new figure so we don't have to . create a dataset just
        % for the test image
        
        htemp = figure();
        newaxes = axes('DataAspectRatio',[1 1 1]);
        colormap(gray);
        subimage(I1,gray(65536));
        %TecPIV_Display(myhandles.TecPivFolder,I1,hPlotAxes,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
        
        hpanelContrast.Visible='off';
        hSecondFigure.Visible = 'off';
         
        guidata(hMainFigure,myhandles);
end
function hPIVsettingsMenuitemCallback (hPIVsettingsMenuitem, eventdata)
        
        % populate the list of sources
        hpopupPIVDataSelector.String = hpopupSourceSelector.String;
        
        % make second figure visible
        hSecondFigure.Visible = 'on';
        hpanelPIVSettings.Visible='on';
        guidata(hMainFigure,myhandles);
end
function hCloseMenuitemCallback (hCloseMenuitem, eventdata)
        disp('Saving figures...')
        saveas(hMainFigure, 'MainFigure.m', 'm');
        saveas(hSecondFigure, 'SecondFigure.m', 'm');
        
        % create structure of handles
        Data = guidata(hMainFigure); 
               
        save(fullfile(myhandles.PathData, myhandles.ProjectID, myhandles.ProjectIDtag),'Data');
        cd(myhandles.TecPivFolder);
        
        s=isempty(gcp('nocreate'));
        if s == 0
            poolobj = gcp('nocreate');
            delete(poolobj);
        end
        disp('Closing now... bye bye')
        diary off
        pause(1);
        fclose('all');
        pause(1);
        eval(['diary ',logname]);
        copyfile(logname,fullfile(myhandles.PathData, myhandles.ProjectID, [ myhandles.ProjectID '.log'] ));
        fclose('all');
        pause(1);
        delete(findobj(0,'type','figure'));
        
end
function hStartPIVCallback(hStartPIV,eventdata)
        %% -- Main PIV function. Reads the parameters from GUI and 
        
        %myhandles.DataSets;
        ThisDataSetNumber=hpopupPIVDataSelector.Value; %get  the selected dataset
        CurrentFrame=str2num(hImgNumber.String); %get the selected frame
        DatasetFolder = myhandles.DataSets{ThisDataSetNumber,1};
        PathData = myhandles.DataSets{ThisDataSetNumber,2};
        ProjectID = myhandles.DataSets{ThisDataSetNumber,3};
        NumberImages = myhandles.DataSets{ThisDataSetNumber,4};
        ImageWidth = myhandles.DataSets{ThisDataSetNumber,5};%%
        ImageHeight = myhandles.DataSets{ThisDataSetNumber,6};%%%
        TimeInc = myhandles.DataSets{ThisDataSetNumber,7};
        
        %FramePath=fullfile(PathData,ProjectID,DatasetFolder,['IMG_' hImgNumber.String '.tif']);
        %I0=imread(FramePath);
        %TecPIV_Display(myhandles.TecPivFolder,I0,hPlotAxes,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
        % Display the first image
        TecPIV_Display_2(hPlotAxes,myhandles.DataSets{ThisDataSetNumber,9},myhandles.DataSets,ThisDataSetNumber,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
        
        % get the parameters from GUI
        
        DataSetNumber=hpopupPIVDataSelector.Value; %get  the selected dataset
        
        StartNumber=str2num(hPIVImgStart.String);
        EndNumber=str2num(hPIVImgEnd.String);
        ImageInc=str2num(hPIVImgInc.String);
        AdjustContrast=hAdjustContrastPIVRadioButton.Value;
        InverseImage = hInverseImageRadioButton.Value;
        SubtractBackgound = hSubtractBackgroundRadioButton.Value;
        GaussianWindowSize = str2num(hGaussianWindowSizeText.String);
        GaussianSigma = str2num(hGaussianSigmaText.String);
        
        UseROI=hUseROIPIVRadioButton.Value;
        DeformROI=hDeformROIPIVRadioButton.Value;
        

        % get processing parameters

        Content=hpopupIntAreaPass1.String;
        IntAreaPass1=str2num(Content{hpopupIntAreaPass1.Value});
        
%         Content=hpopupStepPass1.String;
%         StepPass1=str2num(Content{hpopupStepPass1.Value});
        StepPass1 = IntAreaPass1 /2;
        
        Content=hpopupNumberPass.String;
        NumberPass=str2num(Content{hpopupNumberPass.Value});
        
        Content=hpopupIntAreaPass2.String;
        IntAreaPass2=str2num(Content{hpopupIntAreaPass2.Value});
        
        Content=hpopupIntAreaPass3.String;
        IntAreaPass3=str2num(Content{hpopupIntAreaPass3.Value});
        
        Content=hpopupIntAreaPass4.String;
        IntAreaPass4=str2num(Content{hpopupIntAreaPass4.Value});
        
%         Content=hpopupStepPass2.String;
%         Step2=str2num(Content{hpopupStepPass2.Value});
%         Revert to pivLab way of using 50% overlap for multipass
        Step2 = IntAreaPass2 /2;
        Step3 = IntAreaPass3 /2;
        Step4 = IntAreaPass4 /2;
        
        Content=hpopupWindowDeform.String;
        WindowDeform =  Content{hpopupWindowDeform.Value};
        
        
%         Content=hpopupStepPass3.String;
%         Step3=str2num(Content{hpopupStepPass3.Value});
%         
%         Content=hpopupStepPass4.String;
%         Step4=str2num(Content{hpopupStepPass4.Value});
        
        SubPix=hpopupSubPix.Value;
        
        % get postprocessing values from GUI
        DoVeloFilter=hRadioVelFilter.Value;
        
        umin=str2num(hPIVPostUmin.String);
        umax=str2num(hPIVPostUmax.String);
        vmin=str2num(hPIVPostVmin.String);
        vmax=str2num(hPIVPostVmax.String);
        
        DoStdFilter=hRadioStdFilter.Value;
        ThresStd=str2num(hPIVPostThresStd.String);
        
        DoUnivFilter=hRadioUOFilter.Value;
        EpsMed=str2num(hPIVPostEpsMed.String);
        ThresMed=str2num(hPIVPostThresMed.String);
        
        
        Content=hpopupKernelUO.String;
        KernelUO=Content{hpopupKernelUO.Value};
        
        switch KernelUO
            case '3x3'
                KernelUO=1;
            case '5x5'
                KernelUO=2;
            case '7x7'
                KernelUO=3;
        end
        
        DoInterpolate=hRadioInterpolate.Value;
        Content=hpopupInterpol.String;
        InterpolMethod=Content{hpopupInterpol.Value};
        
        
        param = cell(100,1); 
        
        param{1,1} = IntAreaPass1;
        param{2,1} = StepPass1;
        param{3,1} = SubPix;
        
        param{4,1} = AdjustContrast;
        param{5,1} = InverseImage;
        param{6,1} = SubtractBackgound;
        param{7,1} = GaussianWindowSize;
        param{8,1} = GaussianSigma;
        param{9,1} = UseROI;
        param{10,1} = DeformROI;
        param{11,1} = NumberPass;
        param{12,1} = IntAreaPass2;
        param{13,1} = IntAreaPass3;
        param{14,1} = IntAreaPass4;
        param{15,1} = WindowDeform;
        param{16,1} = umin;
        param{17,1} = umax;
        param{18,1} = vmin;
        param{19,1} = vmax;
        param{20,1} = ThresStd;
        param{21,1} = EpsMed;
        param{22,1} = ThresMed;
        
        param{23,1} = StartNumber;
        param{24,1} = EndNumber;
        param{25,1} = ImageInc;
        
        param{27,1} = Step2;
        param{28,1} = Step3;
        param{29,1} = Step4;
        
        param{30,1} = DoVeloFilter;
        param{31,1} = DoStdFilter;
        
        param{32,1} = DoUnivFilter;
        param{33,1} = KernelUO;
        
        param{34,1} = DoInterpolate;
        param{35,1} = InterpolMethod;
        
        if UseROI == 1
            
            CheckExist = isfield(myhandles ,'ROI');
            if CheckExist == 1 % if mask exist
                % ask user if re-use existing mask
                % Construct a questdlg with 2 options
                ReUseMask = questdlg('A mask already exist. Would you like to re-use it?', ...
                'Yes','No');
                % Handle response
                switch ReUseMask
                    case 'Yes'
                    case 'No'
                        uiwait(msgbox('Using the mouse, specify the region by selecting vertices of the polygon. You can move or resize the polygon using the mouse. Create the mask by double-clicking, or by right-clicking inside the region and selecting Create mask from the context menu.'))
                        [ROI_X, ROI_Y]= ginput(2);
                        ROI_xmin=floor(min(ROI_X));
                        ROI_xmax=ceil(max(ROI_X));
                        ROI_ymin=floor(min(ROI_Y));
                        ROI_ymax=ceil(max(ROI_Y));
                        ROI_width = ROI_xmax - ROI_xmin;
                        ROI_height = ROI_ymax - ROI_ymin;
                        ROI=[ROI_xmin,ROI_ymin,ROI_width,ROI_height];
                        myhandles.ROI=ROI;
                end 
            else % ROI does not exist
                uiwait(msgbox('Using the mouse, specify the region by selecting vertices of the polygon. You can move or resize the polygon using the mouse. Create the mask by double-clicking, or by right-clicking inside the region and selecting Create mask from the context menu.'))
                [ROI_X, ROI_Y]= ginput(2);
                 ROI_xmin=min(ROI_X);
                 ROI_xmax=max(ROI_X);
                 ROI_ymin=min(ROI_Y);
                 ROI_ymax=max(ROI_Y);
                 ROI_width = ROI_xmax - ROI_xmin;
                 ROI_height = ROI_ymax - ROI_ymin;
                 ROI=[ROI_xmin,ROI_ymin,ROI_width,ROI_height];
                 myhandles.ROI=ROI;
            end
            
            param{26,1}=myhandles.ROI;
        end
         
        
        cd( DatasetFolder ); % go into Cam_1/Raw/Rectified
        mkdir( 'Vectors' );  % make folder Cam_1/Raw/Rectified/Vector
        NewDataFolder=fullfile(DatasetFolder,'Vectors'); % create path to new data
        
        cd(fullfile(myhandles.PathData,myhandles.ProjectID)); %go back to main project folder
          
        
        [x,y,u,v,typevector] = TecPIV_Call_PIV(myhandles.DataSets,DataSetNumber,param); % pass all the paremeters assembled in tables to PIV function
        ImageFolderNumber=DataSetNumber;
        %PIVTime=toc;
        %message=sprintf('PIV correlation done in %0.2f s',PIVTime); 
        %disp(message)
        
       
        % Display first image with vector field
        FramePath=fullfile(PathData,ProjectID,DatasetFolder,['IMG_', num2str(StartNumber+ImageInc) ,'.tif']);
        I0=imread(FramePath);
        
        X = [];
        Y = [];
        U = [];
        V = [];
        load(fullfile(PathData,ProjectID,DatasetFolder, 'Vectors',['Vector_' num2str(StartNumber+ImageInc) '.mat']));
        
        myhandles.VectorField{1,1} = 1; % display yes/no
        myhandles.VectorField{1,5} = X; %cell2mat(X); 
        myhandles.VectorField{1,6} = Y; %cell2mat(Y); 
        myhandles.VectorField{1,7} = U; %cell2mat(U); 
        myhandles.VectorField{1,8} = V; %cell2mat(V); 
        myhandles.VectorField{1,9} = [myhandles.PhysU,'/s']; % unit
        myhandles.VectorField{1,10} = TimeInc; % time between successive images in data source
        myhandles.VectorField{1,11} = ImageInc; % increment correlated images
        myhandles.VectorField{1,12} = myhandles.ImScale;
        myhandles.VectorField{1,13} = fullfile(myhandles.ProjectID,NewDataFolder);
        myhandles.VectorField{1,14} = 'black'; % color of vectors 
        myhandles.VectorField{1,15} = typevector;
        myhandles.VectorField{1,16} = 0; % plot as a grid
        
        myhandles.Derivative{1,1} = 1; % display 1=yes 0=no
        
        cla(hPlotAxes);
        %TecPIV_Display(myhandles.TecPivFolder,I0,hPlotAxes,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
        
        %% Add entry to the source selector
        ThisDataSetNumber=myhandles.NumberOfDatasets+1;
        myhandles.NumberOfDatasets=ThisDataSetNumber;

        myhandles.DataSets{ThisDataSetNumber,1} = NewDataFolder;
        myhandles.DataSets{ThisDataSetNumber,2} = PathData;
        myhandles.DataSets{ThisDataSetNumber,3} = myhandles.ProjectID;
        myhandles.DataSets{ThisDataSetNumber,4} = NumberImages-ImageInc;
        myhandles.DataSets{ThisDataSetNumber,5} = [];
        myhandles.DataSets{ThisDataSetNumber,6} = [];
        myhandles.DataSets{ThisDataSetNumber,7} = myhandles.DataSets{DataSetNumber,7}; % same time inc as originial dataset (for now)
        myhandles.DataSets{ThisDataSetNumber,8} = ImageInc;
        myhandles.DataSets{ThisDataSetNumber,9} = StartNumber+ImageInc;
        myhandles.DataSets{ThisDataSetNumber,10} = EndNumber;
        myhandles.DataSets{ThisDataSetNumber,11} = ImageFolderNumber; % the number of the dataset that was used for the PIV calculation and to be plotted under the vector field;
        myhandles.DataSets{ThisDataSetNumber,12} = [];
        myhandles.DataSets{ThisDataSetNumber,13} = [];
        myhandles.DataSets{ThisDataSetNumber,14} = [];
        myhandles.DataSets{ThisDataSetNumber,15} = 'vector';
        myhandles.DataSets{ThisDataSetNumber,16} = 1;
        myhandles.DataSets{ThisDataSetNumber,17} = ThisDataSetNumber;
        
        % Display the first image
        TecPIV_Display_2(hPlotAxes,myhandles.DataSets{ThisDataSetNumber,9},myhandles.DataSets,ThisDataSetNumber,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
        
        NewEntryName=myhandles.DataSets{ThisDataSetNumber,1};
        entries= hpopupSourceSelector.String;
        entries = [entries; NewEntryName];
        hpopupSourceSelector.String = entries;
        hpopupSourceSelector.Value=ThisDataSetNumber; % Select new Entry
        
        
        %% update slider and img number       
         numSteps=floor((EndNumber-StartNumber)/ImageInc);
         hslider.Max=numSteps;
         hslider.Min=1;
         hslider.Value=1; % default shows first vector field
         hslider.SliderStep = [ 1/(numSteps -1),  1/(numSteps-1)]; 
         
         hImgNumber.String=num2str(StartNumber+ImageInc); 
         
         Frame=StartNumber+(hslider.Value-1)*ImageInc;
         CurrentTime = (Frame-StartNumber)*TimeInc; % calculate the time
         hImgTime.String=num2str(CurrentTime); % and time

        % make second figure not visible
        hpanelPIVSettings.Visible='off';
        hSecondFigure.Visible = 'off';
        
        set(hMainFigure, 'pointer','arrow')
        myhandles.DataSets
        
        SavePIVparam=param;
        save('PIV_parameters.mat','SavePIVparam'); 
        
        guidata(hMainFigure,myhandles);
end
function hExportAllFrameMenuitemCallback(hExportAllFrameMenuitem,eventdata)
        % make second figure visible
        hSecondFigure.Visible = 'on';
        hpanelExportAllSettings.Visible='on';
        guidata(hMainFigure,myhandles);
end
function hStartExportAll(hStartExportAllButton,eventdata)
        % start when export button is pressed
        
        % get choices
        DoPDF=hExportAllPDFRadioButton.Value;
        DoPNG=hExportAllPNGRadioButton.Value;
        DoAVI=hExportAllAVIRadioButton.Value;
        
        % get dataset
        SourceNum=hpopupSourceSelector.Value;
        SourceName=char(hpopupSourceSelector.String(SourceNum));

        ThisDataSetNumber=hpopupSourceSelector.Value; % get  the selected dataset
        CurrentFrame=str2num(hImgNumber.String); % get the selected frame
        DatasetFolder = myhandles.DataSets{ThisDataSetNumber,1};
        PathData = myhandles.DataSets{ThisDataSetNumber,2};
        ProjectID = myhandles.DataSets{ThisDataSetNumber,3};
        NumberImages = myhandles.DataSets{ThisDataSetNumber,4};
        ImageWidth = myhandles.DataSets{ThisDataSetNumber,5};
        ImageHeight = myhandles.DataSets{ThisDataSetNumber,6};
        TimeInc = myhandles.DataSets{ThisDataSetNumber,7};
        ImageInc = myhandles.DataSets{ThisDataSetNumber,8};
        StartNumber = myhandles.DataSets{ThisDataSetNumber,9}; 
        EndNumber = myhandles.DataSets{ThisDataSetNumber,10};
       
        % initialise
        X = [];
        Y = [];
        U = [];
        V = [];
        typevector=[];
        
         if DoAVI == 1 % export a movie
             writerObj =  VideoWriter([ProjectID '.avi'],'Uncompressed AVI');
             writerObj.FrameRate=5; % TO DO: Provide way to change that in display settings ?
             open(writerObj);
         end
        
        % check if vector dataset
        k=strfind(DatasetFolder,'Vector');
        
        if isempty(k) == 1 % Dataset is not vector
            
            % Display progress
            progressStepSize = 1;
            ppm = ParforProgMon('Exporting frames: ', NumberImages-1, progressStepSize, 400, 65);
                
            for i=1:NumberImages-1
                Frame=StartNumber+(i)*ImageInc;
                CurrentTime = Frame*TimeInc;
                Framepath=fullfile(PathData,ProjectID,DatasetFolder,['IMG_' num2str(Frame) '.tif']);
                I0=imread(Framepath);
                
                % create new figure
                h = figure();
                newaxes = axes('DataAspectRatio',[1 1 1]);
                
                TecPIV_Display(myhandles.TecPivFolder,I0,newaxes,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
                set(h,'Visible', 'off'); % don't show image when exporting to accelerate
        
                % Set up the paper size / position
                set(h, 'PaperUnits', 'centimeters');
                set(h, 'PaperSize', [21 29.7]);    
                set(h, 'PaperPositionMode', 'manual');
                set(h, 'PaperPosition', [0 0 21 29.7]);
                
                switch SourceName
                    case 'Raw'
                        rootfilename=strcat('Raw-IMG',num2str(Frame));  
                    case 'Calibration'
                        rootfilename=strcat('Calibration-IMG',num2str(Frame)); 
                    case 'Calibration\Rectified'
                        rootfilename=strcat('Calibration-Rectified-IMG',num2str(Frame));
                    case 'Raw\Rectified'
                        rootfilename=strcat('Raw-Rectified-IMG',num2str(Frame));
                    case 'Raw\Rectified\Vectors'
                        rootfilename=strcat('Raw-Rectified-Vectors-IMG',num2str(Frame));
                end
                
                % give title with file name and time
                FrameTitle=strcat(rootfilename,'-',num2str(CurrentTime),'s');
                title(FrameTitle);
                
                if DoPDF == 1 % export PDF (high quality fig)
                    filename=strcat(rootfilename,'.pdf');
                    print(h,'-painters','-dpdf',filename)
                end
                
                if DoPNG == 1 % export png (low quality fig)
                    set(h, 'Units','pixels', 'Position', [0 0 1080 720]);
                    filename=strcat(rootfilename,'.png');
                    print(h,filename,'-dpng')
                end
                
                if DoAVI == 1 % export movie
                    set(h, 'Units','pixels', 'Position', [0 0 1080 720]);
                    F(i)=getframe(h);
                    writeVideo(writerObj,F(i));
                end
                close(h)
                
                % Display progress
                if mod(i,progressStepSize)==0
                     ppm.increment();
                end 

            end
        else %dataset is vector
            
            % check if dataset includes Rectified
            k=strfind(DatasetFolder,'Rectified');
            
            if isempty(k) == 1 % Dataset includes Vector but not Rectified
                ImageFolder=fullfile('Raw');
                % Display progress
                progressStepSize = 1;
                ppm = ParforProgMon('Exporting frames: ', NumberImages-1, progressStepSize, 400, 65);
                for i=1:NumberImages-1
                    Frame=StartNumber+(i-1)*ImageInc;
                    CurrentTime = Frame*TimeInc;
                    Framepath=fullfile(PathData,ProjectID,ImageFolder,['IMG_' num2str(Frame) '.tif']);
                    I0=imread(Framepath);
                    Vector=fullfile(PathData,ProjectID,DatasetFolder,['Vector_' num2str(Frame) '.mat']);
                    load(Vector);
                
                    myhandles.VectorField{1,5} = X; %cell2mat(X); 
                    myhandles.VectorField{1,6} = Y; %cell2mat(Y); 
                    myhandles.VectorField{1,7} = U; %cell2mat(U); 
                    myhandles.VectorField{1,8} = V; %cell2mat(V); 
                
                    % create new figure
                    h = figure();
                    newaxes = axes('DataAspectRatio',[1 1 1]);
                
                    TecPIV_Display(myhandles.TecPivFolder,I0,newaxes,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
                    set(h,'Visible', 'off'); 
        
                    % Set up the paper size / position
                    set(h, 'PaperUnits', 'centimeters');
                    set(h, 'PaperSize', [21 29.7]);    % Set to final desired size here as well as 2 lines below
                    set(h, 'PaperPositionMode', 'manual');
                    set(h, 'PaperPosition', [0 0 21 29.7]);
                    
                    rootfilename = strcat('Unknown-IMG',num2str(Frame));
                
                    switch SourceName
                        case 'Raw'
                            rootfilename=strcat('Raw-IMG',num2str(Frame));  
                        case 'Calibration'
                            rootfilename=strcat('Calibration-IMG',num2str(Frame)); 
                        case 'Calibration\Rectified'
                            rootfilename=strcat('Calibration-Rectified-IMG',num2str(Frame));
                        case 'Raw\Rectified'
                            rootfilename=strcat('Raw-Rectified-IMG',num2str(Frame));
                        case 'Raw\Rectified\Vectors'
                            rootfilename=strcat('Raw-Rectified-Vectors-IMG',num2str(Frame));
                        case 'Raw\Vectors'
                            rootfilename=strcat('Raw-Vectors-IMG',num2str(Frame));
                    end
                    
                FrameTitle=strcat(rootfilename,'-',num2str(CurrentTime),'s');
                title(FrameTitle);
                
                if DoPDF == 1
                    filename=strcat(rootfilename,'.pdf');
                    print(h,'-painters','-dpdf',filename)
                end
                if DoPNG == 1
                    set(h, 'Units','pixels', 'Position', [0 0 1080 720]);
                    filename=strcat(rootfilename,'.png');
                    print(h,filename,'-dpng')
                end
                if DoAVI == 1
                    set(h, 'Units','pixels', 'Position', [0 0 1080 720]);
                    F(i)=getframe(h);
                    writeVideo(writerObj,F(i));
                end
                close(h)
                
                % Display progress
                if mod(i,progressStepSize)==0
                     ppm.increment();
                end 

                end
                
                
                
                
            else  
                % dataset includes vector and Rectified
                ImageFolder=fullfile('Raw','Rectified');
                
                % Display progress
                progressStepSize = 1;
                ppm = ParforProgMon('Exporting frames: ', NumberImages-1, progressStepSize, 400, 65);
                
                for i=1:NumberImages-1
                    Frame=StartNumber+(i-1)*ImageInc;
                    CurrentTime = Frame*TimeInc;
                    Framepath=fullfile(PathData,ProjectID,ImageFolder,['IMG_' num2str(Frame) '.tif']);
                    I0=imread(Framepath);
                    Vector=fullfile(PathData,ProjectID,DatasetFolder,['Vector_' num2str(Frame) '.mat']);
                    load(Vector);
                
                    myhandles.VectorField{1,5} = X; %cell2mat(X); 
                    myhandles.VectorField{1,6} = Y; %cell2mat(Y); 
                    myhandles.VectorField{1,7} = U; %cell2mat(U); 
                    myhandles.VectorField{1,8} = V; %cell2mat(V); 
                
                    % create new figure
                    h = figure();
                    newaxes = axes('DataAspectRatio',[1 1 1]);
                
                    TecPIV_Display(myhandles.TecPivFolder,I0,newaxes,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
                    set(h,'Visible', 'off'); 
        
                    % Set up the paper size / position
                    set(h, 'PaperUnits', 'centimeters');
                    set(h, 'PaperSize', [21 29.7]);    % Set to final desired size here as well as 2 lines below
                    set(h, 'PaperPositionMode', 'manual');
                    set(h, 'PaperPosition', [0 0 21 29.7]);
                    
                    rootfilename=strcat('Unknown-IMG',num2str(Frame));
                
                    switch SourceName
                        case 'Raw'
                            rootfilename=strcat('Raw-IMG',num2str(Frame));  
                        case 'Calibration'
                            rootfilename=strcat('Calibration-IMG',num2str(Frame)); 
                        case 'Calibration\Rectified'
                            rootfilename=strcat('Calibration-Rectified-IMG',num2str(Frame));
                        case 'Raw\Rectified'
                            rootfilename=strcat('Raw-Rectified-IMG',num2str(Frame));
                        case 'Raw\Rectified\Vectors'
                            rootfilename=strcat('Raw-Rectified-Vectors-IMG',num2str(Frame));
                    end
                    
                FrameTitle=strcat(rootfilename,'-',num2str(CurrentTime),'s');
                title(FrameTitle);
                
                if DoPDF == 1
                    filename=strcat(rootfilename,'.pdf');
                    print(h,'-painters','-dpdf',filename)
                end
                if DoPNG == 1
                    set(h, 'Units','pixels', 'Position', [0 0 1080 720]);
                    filename=strcat(rootfilename,'.png');
                    print(h,filename,'-dpng')
                end
                if DoAVI == 1
                    set(h, 'Units','pixels', 'Position', [0 0 1080 720]);
                    F(i)=getframe(h);
                    writeVideo(writerObj,F(i));
                end
                close(h)
                
                % Display progress
                if mod(i,progressStepSize)==0
                     ppm.increment();
                end 

                end
                
                if DoAVI == 1
                close(writerObj);
                end      
            end    
        end
            
           % make second figure visible
           hpanelExportAllSettings.Visible='off'; 
        hSecondFigure.Visible = 'off';
        guidata(hMainFigure,myhandles);
end

% function hPlayCallback(hplay,eventdata)
%         % Function returned when the play button is pushed. Start the
%         % animation in axes
%         
%         ThisDataSetNumber=hpopupSourceSelector.Value; %get  the selected dataset
%         CurrentFrame=str2num(hImgNumber.String); %get the selected frame
%         DatasetFolder = myhandles.DataSets{ThisDataSetNumber,1};
%         PathData = myhandles.DataSets{ThisDataSetNumber,2};
%         ProjectID = myhandles.DataSets{ThisDataSetNumber,3};
%         NumberImages = myhandles.DataSets{ThisDataSetNumber,4};
%         ImageWidth = myhandles.DataSets{ThisDataSetNumber,5};%%
%         ImageHeight = myhandles.DataSets{ThisDataSetNumber,6};%%%
%         TimeInc = myhandles.DataSets{ThisDataSetNumber,7};
%         ImageInc = myhandles.DataSets{ThisDataSetNumber,8};
%         StartNumber = myhandles.DataSets{ThisDataSetNumber,9}; % This is because the start number may have been changed when doing the correlation
%         EndNumber = myhandles.DataSets{ThisDataSetNumber,10}; % Same here, may have been adjusted during correlation procedure
%         
%         myhandles.Stop= 0;
%         Frame = CurrentFrame;
%         
%          % check if vector dataset
%         k=strfind(DatasetFolder,'Vector');
%        
%         if isempty(k) == 1 % Dataset is not vector
%             myhandles.VectorField{1,1} = 0; % display yes/no
%             myhandles.Derivative{1,1} = 0;
%             
%             while myhandles.Stop == 0
%                 
%                 CurrentTime = Frame*TimeInc;
%                 Framepath=fullfile(PathData,ProjectID,ImageFolder,['IMG_' num2str(Frame) '.tif']);
%                 
% 
%                 cla(hPlotAxes);
%                 TecPIV_Display(myhandles.TecPivFolder,I0,hPlotAxes,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative)
%                 drawnow
%                 
%                 
%                 Frame=Frame+ImageInc;
%                 
%                 if Frame > EndNumber
%                     Frame = StartNumber;
%                 end
%                  
%             end
%             
%         else % Dataset is vector
%              myhandles.VectorField{1,1} = 1; % display vectors
%              
%             % check if dataset includes Rectified
%             k=strfind(DatasetFolder,'Rectified');
%             
%             if isempty(k) == 1 % Dataset includes Vector but not Rectified
%                 ImageFolder=fullfile('Raw');
%             else  
%                 % dataset includes vector and Rectified
%                 ImageFolder=fullfile('Raw','Rectified');   
%             end
%             
%             while myhandles.Stop == 0
%                 
%                 Vector=fullfile(PathData,ProjectID,DatasetFolder,['Vector_' num2str(Frame) '.mat']);
%                 % initialise
%                 X = [];
%                 Y = [];
%                 U = [];
%                 V = [];
%                 
%                 load(Vector);
%                 myhandles.VectorField{1,5} = X; 
%                 myhandles.VectorField{1,6} = Y; 
%                 myhandles.VectorField{1,7} = U;  
%                 myhandles.VectorField{1,8} = V; 
%                 
%                
%                 CurrentTime = Frame*TimeInc;
%                 
%                 hImgNumber.String=num2str(Frame); % update img num
%                 hImgTime.String=num2str(CurrentTime); % and time
%                 
%                 hslider.Value=1+(Frame-StartNumber)/ImageInc;
%                 
%                 Framepath=fullfile(PathData,ProjectID,ImageFolder,['IMG_' num2str(Frame) '.tif']);
%                 I0=imread(Framepath);
%                 cla(hPlotAxes);
%                 
%                 TecPIV_Display(myhandles.TecPivFolder,I0,hPlotAxes,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative)
%                 drawnow
% 
%                 Frame=Frame+ImageInc;
%                 
%                 if Frame > EndNumber
%                     Frame= StartNumber;
%                 end
%                  
%             end
%             
%                 
%             
%         end
%         
%         guidata(hMainFigure,myhandles);
% end

%function hStopCallback(hstop,eventdata)
        %myhandles.Stop = 1;
%end

function hSourceSelectorCallback(hpopupSourceSelector,eventdata)
        % called when selected value in source selector is changed
        
        % get the selected value & name
        SourceNum=hpopupSourceSelector.Value;
        SourceName=char(hpopupSourceSelector.String(SourceNum));
        
        ThisDataSetNumber=SourceNum; % the selected dataset
        CurrentFrame=str2num(hImgNumber.String); % the selected frame
        
        DatasetFolder = myhandles.DataSets{ThisDataSetNumber,1}; % dataset name
        PathData = myhandles.DataSets{ThisDataSetNumber,2}; % TecPIV folder
        ProjectID = myhandles.DataSets{ThisDataSetNumber,3}; % project name
        NumberImages = myhandles.DataSets{ThisDataSetNumber,4}; %number images
        ImageWidth = myhandles.DataSets{ThisDataSetNumber,5};% image width
        ImageHeight = myhandles.DataSets{ThisDataSetNumber,6};% image height
        TimeInc = myhandles.DataSets{ThisDataSetNumber,7}; % time inc
        ImageInc = myhandles.DataSets{ThisDataSetNumber,8}; % image inc    
        StartImage = myhandles.DataSets{ThisDataSetNumber,9}; % StartImage
        EndImage = myhandles.DataSets{ThisDataSetNumber,10}; % EndImage
        
        %FramePath=fullfile(PathData,ProjectID,DatasetFolder,['IMG_' num2str(StartImage) '.tif']);
        %I0=imread(FramePath);
        %TecPIV_Display(myhandles.TecPivFolder,I0,hPlotAxes,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
        
        % Display the first image
        TecPIV_Display_2(hPlotAxes,StartImage,myhandles.DataSets,ThisDataSetNumber,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
        
        
        % update the player
        [ myhandles.ImgNumber, myhandles.Time, myhandles.SliderMin, myhandles.SliderMax,myhandles.SliderStep , myhandles.SliderVal] = TecPIV_UpdatePlayer('FrameNum',...
            num2str(StartImage),...
            myhandles.DataSets,...
            ThisDataSetNumber);
         
        if NumberImages == 1
            hslider.Visible='off';
            hImgNumber.Visible='off';
            hImgTime.Visible='off';
        else
            hImgNumber.String=num2str(myhandles.ImgNumber);
            hImgTime.String=num2str(myhandles.Time);
            hslider.SliderStep=myhandles.SliderStep;
            hslider.Max=myhandles.SliderMax;
            hslider.Min=myhandles.SliderMin;
            hslider.Value=str2double(myhandles.SliderVal);
        end
        
        guidata(hMainFigure,myhandles); 
              
        
end   
function hSliderCallback(hslider,eventdata)
       % This function is activated when the slider has been triggered 
       % Update the image number, image and time...
             
       ThisDataSetNumber=hpopupSourceSelector.Value; % get  the selected dataset
       CurrentFrame=str2num(hImgNumber.String); % get the selected frame
       DatasetFolder = myhandles.DataSets{ThisDataSetNumber,1};
       PathData = myhandles.DataSets{ThisDataSetNumber,2};
       ProjectID = myhandles.DataSets{ThisDataSetNumber,3};
       NumberImages = myhandles.DataSets{ThisDataSetNumber,4};
       ImageWidth = myhandles.DataSets{ThisDataSetNumber,5};%%
       ImageHeight = myhandles.DataSets{ThisDataSetNumber,6};%%%
       TimeInc = myhandles.DataSets{ThisDataSetNumber,7};
       ImageInc = myhandles.DataSets{ThisDataSetNumber,8};
       StartNumber = myhandles.DataSets{ThisDataSetNumber,9}; % This is because the start number may have been changed when doing the correlation
       EndNumber = myhandles.DataSets{ThisDataSetNumber,10}; % Same here, may have been adjusted during correlation procedure
       
       NewSliderValue=round(hslider.Value);
       
       if ceil(hslider.Value) == floor(hslider.Value) % is integer
           NewSliderValue=hslider.Value;
       else
           NewSliderValue=round(hslider.Value);
       end
       NewFrame=StartNumber+(NewSliderValue-1)*ImageInc; % new frame value
        
       CurrentTime = (NewFrame-StartNumber)*TimeInc; % calculate the time
       hImgNumber.String=num2str(NewFrame); % update img num
       hImgTime.String=num2str(CurrentTime); % and time
       
       % Now update the image
       
        k=strfind(DatasetFolder,'Vector'); % check if vector dataset
       
        if isempty(k) == 1 % Dataset is not vector
            myhandles.VectorField{1,1} = 0; % display vectors 0=yes 1=no
            myhandles.Derivative{1,1} = 0; % display derivatives 0=yes 1=no
            Framepath=fullfile(PathData,ProjectID,DatasetFolder,['IMG_' num2str(NewFrame) '.tif']);
            I0=imread(Framepath);
            cla(hPlotAxes);
            %TecPIV_Display(myhandles.TecPivFolder,I0,hPlotAxes,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative)
            % Display the first image
            TecPIV_Display_2(hPlotAxes,NewFrame,myhandles.DataSets,ThisDataSetNumber,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
            drawnow
            
        else % dataset is vector 
            k=strfind(DatasetFolder,'Rectified'); % check if dataset includes Rectified
            
            if isempty(k) == 1 % Dataset is vector, not Rectified
            ImageFolder=fullfile('Raw');          
            
            else % Dataset is vector, includes Rectified      
            ImageFolder=fullfile('Raw','Rectified');
            end
                    
            Framepath=fullfile(PathData,ProjectID,ImageFolder,['IMG_' num2str(NewFrame) '.tif']);
            I0=imread(Framepath);
            
            Vector=fullfile(PathData,ProjectID,DatasetFolder,['Vector_' num2str(NewFrame) '.mat']);
            Data=load(Vector);
            X=Data.X;
            Y=Data.Y;
            U=Data.U;
            V=Data.V;
            
            myhandles.VectorField{1,1} = 1; % display yes/no
            myhandles.VectorField{1,5} = X; %cell2mat(X); 
            myhandles.VectorField{1,6} = Y; %cell2mat(Y); 
            myhandles.VectorField{1,7} = U; %cell2mat(U); 
            myhandles.VectorField{1,8} = V; %cell2mat(V); 
            
            myhandles.VectorField{1,9} = [myhandles.PhysU,'/s']; % unit
            myhandles.VectorField{1,10} = TimeInc; % time between successive images in data source
            myhandles.VectorField{1,11} = ImageInc; % increment correlated images
            myhandles.VectorField{1,12} = myhandles.ImScale;
            myhandles.VectorField{1,13} = fullfile(myhandles.ProjectID,DatasetFolder);
            myhandles.VectorField{1,14} = 'black'; % color of vectors
            %myhandles.VectorField{1,15} = typevector;

            
            cla(hPlotAxes);
            %TecPIV_Display(myhandles.TecPivFolder,I0,hPlotAxes,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative)
            % Display the first image
            TecPIV_Display_2(hPlotAxes,NewFrame,myhandles.DataSets,ThisDataSetNumber,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
            drawnow
            %drawnow
            
            
        end
                

end   
function hImgNumCallback(hImgNum,eventdata)
        % ImgNum has been changed
        Frame=str2num(hImgNumber.String);
        
        ThisDataSetNumber=hpopupSourceSelector.Value; %get  the selected dataset
        CurrentFrame=str2num(hImgNumber.String); %get the selected frame
        DatasetFolder = myhandles.DataSets{ThisDataSetNumber,1};
        PathData = myhandles.DataSets{ThisDataSetNumber,2};
        ProjectID = myhandles.DataSets{ThisDataSetNumber,3};
        NumberImages = myhandles.DataSets{ThisDataSetNumber,4};
        ImageWidth = myhandles.DataSets{ThisDataSetNumber,5};%%
        ImageHeight = myhandles.DataSets{ThisDataSetNumber,6};%%%
        TimeInc = myhandles.DataSets{ThisDataSetNumber,7};
        ImageInc = myhandles.DataSets{ThisDataSetNumber,8};
        StartNumber = myhandles.DataSets{ThisDataSetNumber,9}; % This is because the start number may have been changed when doing the correlation
        EndNumber = myhandles.DataSets{ThisDataSetNumber,10}; % Same here, may have been adjusted during correlation procedure
        
        % check that requested frame is within limits. adjust if necessary
        if Frame > EndNumber
            Frame=EndNumber;
        end
        
        if Frame < StartNumber
            Frame=StartNumber;
        end
        
        % check 
        % if request frame number that is not available, find the closest
        % and adjust
        
        RqFrame=Frame;
        FrameNumbers=[StartNumber:ImageInc:EndNumber]; %list of available frame numbers
        DiffRqNums=abs(FrameNumbers-RqFrame); %difference between rquested and available
        BestNumPos=min(DiffRqNums); % closest available frame
        Frame=FrameNumbers(DiffRqNums==BestNumPos); %set frame number to closest available
        
        % little test to check that we are not at equal distance from 2
        % possible values. If yes then take lower one.
        [n,m]=size(Frame);
        if m==2
            Frame=Frame(1,1);
        end
        
        CurrentTime = (Frame-StartNumber)*TimeInc; % calculate the time
        hImgTime.String=num2str(CurrentTime); % and time
        hImgNum.String=num2str(Frame);
        hslider.Value=1+(Frame-StartNumber)/ImageInc;
        
        % Now update the image
        % check if vector dataset
        k=strfind(DatasetFolder,'Vector');
       
        if isempty(k) == 1 % Dataset is not vector
            myhandles.VectorField{1,1} = 0; % display yes/no
            myhandles.Derivative{1,1} = 0;
            Framepath=fullfile(PathData,ProjectID,DatasetFolder,['IMG_' num2str(Frame) '.tif']);
            I0=imread(Framepath);
            cla(hPlotAxes);
            %TecPIV_Display(myhandles.TecPivFolder,I0,hPlotAxes,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative)
            TecPIV_Display_2(hPlotAxes,Frame,myhandles.DataSets,ThisDataSetNumber,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
            drawnow
            %drawnow
        else % data is vector
            k=strfind(DatasetFolder,'Rectified'); % check if dataset includes Rectified
            
            if isempty(k) == 1 % Dataset is vector, not Rectified
            ImageFolder=fullfile('Raw');          
            
            else % Dataset is vector, includes Rectified      
            ImageFolder=fullfile('Raw','Rectified');
            end
                    
            Framepath=fullfile(PathData,ProjectID,ImageFolder,['IMG_' num2str(Frame) '.tif']);
            I0=imread(Framepath);
            
            Vector=fullfile(PathData,ProjectID,DatasetFolder,['Vector_' num2str(Frame) '.mat']);
            Data=load(Vector);
            X=Data.X;
            Y=Data.Y;
            U=Data.U;
            V=Data.V;
            
            myhandles.VectorField{1,1} = 1; % display yes/no
            myhandles.VectorField{1,5} = X; %cell2mat(X); 
            myhandles.VectorField{1,6} = Y; %cell2mat(Y); 
            myhandles.VectorField{1,7} = U; %cell2mat(U); 
            myhandles.VectorField{1,8} = V; %cell2mat(V); 
            
            myhandles.VectorField{1,9} = [myhandles.PhysU,'/s']; % unit
            myhandles.VectorField{1,10} = TimeInc; % time between successive images in data source
            myhandles.VectorField{1,11} = ImageInc; % increment correlated images
            myhandles.VectorField{1,12} = myhandles.ImScale;
            myhandles.VectorField{1,13} = fullfile(myhandles.ProjectID,DatasetFolder);
            myhandles.VectorField{1,14} = 'black'; % color of vectors 
%            myhandles.VectorField{1,15} = typevector;

            
            cla(hPlotAxes);
            %TecPIV_Display(myhandles.TecPivFolder,I0,hPlotAxes,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative)
            %drawnow
            TecPIV_Display_2(hPlotAxes,Frame,myhandles.DataSets,ThisDataSetNumber,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
            drawnow
            
        end
 
end   
function hImgTimeCallback(hImgTime,eventdata)
        RqTime=str2double(hImgTime.String);
        
        ThisDataSetNumber=hpopupSourceSelector.Value; %get  the selected dataset
        CurrentFrame=str2num(hImgNumber.String); %get the selected frame
        DatasetFolder = myhandles.DataSets{ThisDataSetNumber,1};
        PathData = myhandles.DataSets{ThisDataSetNumber,2};
        ProjectID = myhandles.DataSets{ThisDataSetNumber,3};
        NumberImages = myhandles.DataSets{ThisDataSetNumber,4};
        ImageWidth = myhandles.DataSets{ThisDataSetNumber,5};%%
        ImageHeight = myhandles.DataSets{ThisDataSetNumber,6};%%%
        TimeInc = myhandles.DataSets{ThisDataSetNumber,7};
        ImageInc = myhandles.DataSets{ThisDataSetNumber,8};
        StartNumber = myhandles.DataSets{ThisDataSetNumber,9}; % This is because the start number may have been changed when doing the correlation
        EndNumber = myhandles.DataSets{ThisDataSetNumber,10}; % Same here, may have been adjusted during correlation procedure
        
        FrameNumbers=[StartNumber:ImageInc:EndNumber]; %list of available frame numbers
        RqNum=RqTime/TimeInc+StartNumber; % calculate requested frame number
        
        DiffRqNums=abs(FrameNumbers-RqNum); %difference between rquested and available
        BestNumPos=min(DiffRqNums); % closest available frame
        Frame=FrameNumbers(DiffRqNums==BestNumPos); %set frame number to closest available
        
        % little test to check that we are not at equal distance from 2
        % possible values. If yes then take lower one.
        [n,m]=size(Frame);
        if m==2
            Frame=Frame(1,1);
        end
        CurrentTime = (Frame-StartNumber)*TimeInc; % calculate the time
        hImgTime.String=num2str(CurrentTime); % and time
        hImgNumber.String=num2str(Frame);
        hslider.Value=1+(Frame-StartNumber)/ImageInc;
        
        % Now update the image
        % check if vector dataset
        k=strfind(DatasetFolder,'Vector');
       
        if isempty(k) == 1 % Dataset is not vector
            myhandles.VectorField{1,1} = 0; % display yes/no
            myhandles.Derivative{1,1} = 0;
            Framepath=fullfile(PathData,ProjectID,DatasetFolder,['IMG_' num2str(Frame) '.tif']);
            I0=imread(Framepath);
            cla(hPlotAxes);
            %TecPIV_Display(myhandles.TecPivFolder,I0,hPlotAxes,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative)
            %drawnow
            TecPIV_Display_2(hPlotAxes,Frame,myhandles.DataSets,ThisDataSetNumber,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
            drawnow
        else
            % data is vector
            k=strfind(DatasetFolder,'Rectified'); % check if dataset includes Rectified
            
            if isempty(k) == 1 % Dataset is vector, not Rectified
            ImageFolder=fullfile('Raw');          
            
            else % Dataset is vector, includes Rectified      
            ImageFolder=fullfile('Raw','Rectified');
            end
                    
            Framepath=fullfile(PathData,ProjectID,ImageFolder,['IMG_' num2str(Frame) '.tif']);
            I0=imread(Framepath);
            
            Vector=fullfile(PathData,ProjectID,DatasetFolder,['Vector_' num2str(Frame) '.mat']);
            Data=load(Vector);
            X=Data.X;
            Y=Data.Y;
            U=Data.U;
            V=Data.V;
            
            myhandles.VectorField{1,1} = 1; % display yes/no
            myhandles.VectorField{1,5} = X; %cell2mat(X); 
            myhandles.VectorField{1,6} = Y; %cell2mat(Y); 
            myhandles.VectorField{1,7} = U; %cell2mat(U); 
            myhandles.VectorField{1,8} = V; %cell2mat(V); 
            
            myhandles.VectorField{1,9} = [myhandles.PhysU,'/s']; % unit
            myhandles.VectorField{1,10} = TimeInc; % time between successive images in data source
            myhandles.VectorField{1,11} = ImageInc; % increment correlated images
            myhandles.VectorField{1,12} = myhandles.ImScale;
            myhandles.VectorField{1,13} = fullfile(myhandles.ProjectID,DatasetFolder);
            myhandles.VectorField{1,14} = 'black'; % color of vectors 
            myhandles.VectorField{1,15} = typevector;

            
            cla(hPlotAxes);
            %TecPIV_Display(myhandles.TecPivFolder,I0,hPlotAxes,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative)
            %drawnow
            TecPIV_Display_2(hPlotAxes,Frame,myhandles.DataSets,ThisDataSetNumber,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
            drawnow
        end
  
end
function hApplyDisplaySetting(hApplyDisplaySettingButton,eventdata)
        % read the data in GUI and place in variables
        % background image
        % color palette
        ListCPT=hpopupBackgroundCPTSelector.String;
        myhandles.RawCpt{1,1} = char(ListCPT(hpopupBackgroundCPTSelector.Value));
        myhandles.RawCpt{1,2} = hMinBackgroundColorPalette.String;
        myhandles.RawCpt{1,3} = hMaxBackgroundColorPalette.String;
        
        % vector field
        myhandles.VectorField{1,1} = hDisplayVectorRadioButton.Value;
        myhandles.VectorField{1,2} = [hpopupVecGridFactorSelector.Value, hpopupVecGridFactorSelector.Value];
        
        temp = hpopupVecColorScalingModeSelector.Value;
        if temp == 1
            myhandles.VectorField{1,3} = 'median';
        elseif temp == 2
            myhandles.VectorField{1,3} = 'max';
        else
            myhandles.VectorField{1,3} = 'manual';
        end
 
        myhandles.VectorField{1,4} = str2num(hVecScalingLength.String); % reference value when manual
        
        temp = hpopupVecColorSelector.Value;
        TempList = hpopupVecColorSelector.String;
        myhandles.VectorField{1,14} = TempList(hpopupVecColorSelector.Value); % color of vector field
        
        % derivative
        myhandles.Derivative{1,1} = hDisplayVectorDerivativeRadioButton.Value; % display yes/no
        
        ListPopUpDeriv =  hpopupVecDerivativeTypeSelector.String;
        myhandles.Derivative{1,2} = char(ListPopUpDeriv(hpopupVecDerivativeTypeSelector.Value)); % type of derivative
        
        myhandles.Derivative{1,3} = hpopupVecDerivativeDisplayRangeSelector.Value; %1 = minmax; 2 = +- max; 3 = custom
        
        
        ListCPTDeriv=hpopupVecDerivativeCPTSelector.String; % list of CPT for derivative
        myhandles.Derivative{1,4} = strtrim(char(ListCPTDeriv(hpopupVecDerivativeCPTSelector.Value)));
        
        myhandles.Derivative{1,5} = str2double(hMinVecDerivativeColorPalette.String); 
        myhandles.Derivative{1,6} = str2double(hMaxVecDerivativeColorPalette.String);
        
        myhandles.Derivative{1,7} = str2double(hAlphaDeriv.String); % alpha derivative
        myhandles.Derivative{1,8} = hInterpDerivativeRadioButton.Value; % toggle interp
        ListInterp=hpopupVecDerivativeMethodSelector.String; % list of methods
        myhandles.Derivative{1,9} = char(ListInterp(hpopupVecDerivativeMethodSelector.Value)); %text corresponding to selected method
        
     
        % make second figure and its panels invisible after
        hpanelChangeBackgroundDisplaySettings.Visible='off';
        hpanelChangeVectorDisplaySettings.Visible='off';
        hpanelChangeVectorDerivativeDisplaySettings.Visible='off';
        hpanelChangeApplyDisplaySettings.Visible = 'off';
        hSecondFigure.Visible = 'off';
                
        X = [];
        Y = [];
        U = [];
        V = [];
        typevector=[];
        
        guidata(hMainFigure,myhandles);
        
        % update the figure
        % get figure properties
        ThisDataSetNumber=hpopupSourceSelector.Value; %get  the selected dataset
        CurrentFrame=str2num(hImgNumber.String); %get the selected frame
        DatasetFolder = myhandles.DataSets{ThisDataSetNumber,1};
        PathData = myhandles.DataSets{ThisDataSetNumber,2};
        ProjectID = myhandles.DataSets{ThisDataSetNumber,3};
        
        % check if vector dataset
        k=strfind(DatasetFolder,'Vector');
        
        if isempty(k) == 1 % Dataset is not vector
            Framepath=fullfile(PathData,ProjectID,DatasetFolder,['IMG_' num2str(CurrentFrame) '.tif']);
            I0=imread(Framepath);
        else
            % check if dataset includes Rectified
            k=strfind(DatasetFolder,'Rectified');
            
            if isempty(k) == 1 % Dataset includes Vector but not Rectified
                % find the Cam number
                %[token remain]=strtok(DatasetFolder,'\');
                ImageFolder=fullfile('Raw');
                 Framepath=fullfile(PathData,ProjectID,ImageFolder,['IMG_' num2str(CurrentFrame) '.tif']);
                I0=imread(Framepath);
                Vector=fullfile(PathData,ProjectID,DatasetFolder,['Vector_' num2str(CurrentFrame) '.mat']);
                load(Vector);
                
                myhandles.VectorField{1,5} = X; %cell2mat(X); 
                myhandles.VectorField{1,6} = Y; %cell2mat(Y); 
                myhandles.VectorField{1,7} = U; %cell2mat(U); 
                myhandles.VectorField{1,8} = V; %cell2mat(V);
                myhandles.VectorField{1,15} = typevector; %cell2mat(V);
                
            else  
                % dataset includes vector and Rectified
                % find the Cam number
                %[token remain]=strtok(DatasetFolder,'\');
                ImageFolder=fullfile('Raw','Rectified');
                
                Framepath=fullfile(PathData,ProjectID,ImageFolder,['IMG_' num2str(CurrentFrame) '.tif']);
                I0=imread(Framepath);
                Vector=fullfile(PathData,ProjectID,DatasetFolder,['Vector_' num2str(CurrentFrame) '.mat']);
                load(Vector);
                
                myhandles.VectorField{1,5} = X; %cell2mat(X); 
                myhandles.VectorField{1,6} = Y; %cell2mat(Y); 
                myhandles.VectorField{1,7} = U; %cell2mat(U); 
                myhandles.VectorField{1,8} = V; %cell2mat(V); 
                myhandles.VectorField{1,15} = typevector;
                
                
            end
        end
        
        %TecPIV_Display(myhandles.TecPivFolder,I0,hPlotAxes,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
        TecPIV_Display_2(hPlotAxes,CurrentFrame,myhandles.DataSets,ThisDataSetNumber,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
            drawnow
        
        guidata(hMainFigure,myhandles); 
        
end  
function hLagrangianSumMenuitemCallback(hLagrangianSumMenuitem,eventdata)
        
        ThisDataSetNumber=hpopupSourceSelector.Value; % get  the selected dataset
        SourceData = myhandles.DataSets{ThisDataSetNumber,1};
        NumberSourceData=ThisDataSetNumber;
        
        TecPIV_Lag_Sum(myhandles.DataSets,ThisDataSetNumber)
        
        
        NumberDatasets = myhandles.NumberOfDatasets;
        
        NLag = 0;
        for j = 1:NumberDatasets
            DatasetName =  myhandles.DataSets{j,1};
            TF = contains(DatasetName,'Lagrangian');
            if TF == 1
                NLag = NLag+1;
            end  
        end
        
        % create new dataset
        ThisDataSetNumber=myhandles.NumberOfDatasets+1;
        
        myhandles.NumberOfDatasets=ThisDataSetNumber;
        ThisDataSetName=[myhandles.DataSets{NumberSourceData,1} '/Lagrangian_Sum_' num2str(NLag+1)]; %'Raw/Vectors/Lagrangian_Sum';
        myhandles.DataSets{ThisDataSetNumber,1}=ThisDataSetName;
        myhandles.DataSets{ThisDataSetNumber,2}=myhandles.PathData;
        myhandles.DataSets{ThisDataSetNumber,3}=myhandles.ProjectID;
        myhandles.DataSets{ThisDataSetNumber,4}=myhandles.DataSets{NumberSourceData,4};
        myhandles.DataSets{ThisDataSetNumber,5}=myhandles.DataSets{NumberSourceData,5};
        myhandles.DataSets{ThisDataSetNumber,6}=myhandles.DataSets{NumberSourceData,6};
        myhandles.DataSets{ThisDataSetNumber,7}=myhandles.DataSets{NumberSourceData,7};
        myhandles.DataSets{ThisDataSetNumber,8}=myhandles.DataSets{NumberSourceData,8};
        myhandles.DataSets{ThisDataSetNumber,9}=myhandles.DataSets{NumberSourceData,9};
        myhandles.DataSets{ThisDataSetNumber,10}=myhandles.DataSets{NumberSourceData,10};
        myhandles.DataSets{ThisDataSetNumber,11}=myhandles.DataSets{NumberSourceData,11};
        
        myhandles.DataSets{ThisDataSetNumber,12}=myhandles.DataSets{NumberSourceData,12};
        myhandles.DataSets{ThisDataSetNumber,13}=myhandles.DataSets{NumberSourceData,13};
        myhandles.DataSets{ThisDataSetNumber,14}=myhandles.DataSets{NumberSourceData,14};
        myhandles.DataSets{ThisDataSetNumber,15}='Lagrangian';
        
        myhandles.DataSets{ThisDataSetNumber,16}=myhandles.DataSets{NumberSourceData,16};
        myhandles.DataSets{ThisDataSetNumber,17}=myhandles.DataSets{NumberSourceData,17}; 
        
        myhandles.DataSets{ThisDataSetNumber,18} = NumberSourceData;
        
        myhandles.entries = hpopupSourceSelector.String;
        myhandles.entries = [myhandles.entries; ThisDataSetName];
        hpopupSourceSelector.String = myhandles.entries;
        myhandles.SelectedEntry=ThisDataSetNumber;
        hpopupSourceSelector.Value=ThisDataSetNumber; 
        
        Frame=myhandles.DataSets{ThisDataSetNumber,9};
        SaveDatasets=myhandles.DataSets
        save('Datasets.mat','SaveDatasets');
        
        % set derivative to finite strain
        myhandles.Derivative{1,2} = 'Exy';
        
        
        guidata(hMainFigure,myhandles);
    
        % Display the first image
        TecPIV_Display_2(hPlotAxes,Frame,myhandles.DataSets,ThisDataSetNumber,myhandles.RawCpt,myhandles.VectorField,myhandles.Derivative);
        
        
        
%         
%         ImageFolderNumber=myhandles.DataSets{ThisDataSetNumber,11};
%         ImageFolder=myhandles.DataSets{ImageFolderNumber,1};
%         Framepath=fullfile(myhandles.PathData,myhandles.ProjectID,ImageFolder,['IMG_' num2str(Frame) '.tif']);
%         I=imread(Framepath);
% 
%         DatasetFolder = ThisDataSetName;
%         
%         
%         Vector=fullfile(myhandles.PathData,myhandles.ProjectID,DatasetFolder,['Vector_' num2str(Frame) '.mat']);
%         Data=load(Vector);
%         X=Data.X;
%         Y=Data.Y;
%         U=Data.U;
%         V=Data.V;
%         
%         [theta,r] = cart2pol(U,V);
%         Ax = hPlotAxes;
%         axes(Ax);
%         cla(Ax); 
%         subimage(I,gray(65536));
%         h2 = pcolor(X,Y,r);
%         set(h2,'EdgeColor','white');
        
        guidata(hMainFigure,myhandles);
end
function hEulerianSumMenuitemCallback(hEulerianSumMenuitem,eventdata)
    ThisDataSetNumber=hpopupSourceSelector.Value; % get  the selected dataset
    %SourceData = myhandles.DataSets{ThisDataSetNumber,1};
    TecPIV_Cumulative_Eul(myhandles.DataSets,ThisDataSetNumber) 
end
 

guidata(hMainFigure,myhandles); 
end
