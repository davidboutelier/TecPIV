function TracerInit = AddLagrangianTracer(NTracers, TracerInit)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

pad = 5; % padding

% default frame number
DFN = 2;

NextEntry = NTracers + 1; 
hFigAddLagTracer = figure('position',[100 100 800 600],...
    'name', 'Lagrangian tracers');

hApplyTracer = uicontrol('style', 'pushbutton',...
    'parent', hFigAddLagTracer,...
    'string', 'Apply',...
    'position', [(hFigAddLagTracer.Position(3) - 60) 10 50 40]);

h1 = uicontrol('style','text',...
    'parent', hFigAddLagTracer,...
    'string', 'Tracer number :',...
    'HorizontalAlignment', 'left',...
    'position', [10 50 200 30]);

h1.Position(3) = h1.Extent(3) + pad;
h1.Position(4) = h1.Extent(4) + pad;

uiTracerNum = uicontrol('style','edit',...
     'parent', hFigAddLagTracer,...
     'string', num2str(NextEntry),...
     'position', [(h1.Position(1)+h1.Position(3)) 50 30 30]);

uiTracerNum.Position(3) = uiTracerNum.Extent(3) + pad;
uiTracerNum.Position(4) = uiTracerNum.Extent(4) + pad;

h2 = uicontrol('style','text',...
    'parent', hFigAddLagTracer,...
    'string', 'Frame number :',...
    'HorizontalAlignment', 'left',...
    'position', [(uiTracerNum.Position(1)+uiTracerNum.Position(3) + 4 * pad) 50 200 30]);

h2.Position(3) = h2.Extent(3) + pad;
h2.Position(4) = h2.Extent(4) + pad;

uiFrameNum = uicontrol('style','edit',...
     'parent', hFigAddLagTracer,...
     'string', num2str(DFN),...
     'position', [(h2.Position(1)+h2.Position(3)) 50 30 30]);

uiFrameNum.Position(3) = uiFrameNum.Extent(3) + pad;
uiFrameNum.Position(4) = uiFrameNum.Extent(4) + pad;

uibuttoncoord = uicontrol('style', 'pushbutton',...
    'parent', hFigAddLagTracer,...
    'string', 'set coordinates',...
    'position', [(uiFrameNum.Position(1)+uiFrameNum.Position(3) + 4 * pad) 50 200 30]);

uibuttoncoord.Position(3) = uibuttoncoord.Extent(3) + pad;
uibuttoncoord.Position(4) = uibuttoncoord.Extent(4) + pad;

uipopmenuparam = uicontrol('style', 'popupmenu',...
    'parent', hFigAddLagTracer,...
    'string', {'u','v','du/dx','du/dy','vorticity'},...
    'value', 1, ...
    'position', [(uibuttoncoord.Position(1)+uibuttoncoord.Position(3) + 4 * pad) 50 200 30]);

uipopmenuparam.Position(3) = 1.5 * uibuttoncoord.Position(3);
uipopmenuparam.Position(4) = uibuttoncoord.Position(4);

uibuttonparam = uicontrol('style', 'pushbutton',...
    'parent', hFigAddLagTracer,...
    'string', 'add parameter',...
    'position', [(uipopmenuparam.Position(1)+uipopmenuparam.Position(3) + 0 * pad) 50 200 30]);

uibuttonparam.Position(3) = uibuttonparam.Extent(3) + pad;
uibuttonparam.Position(4) = uibuttonparam.Extent(4) + pad;
uibuttonparam.Callback = @AddParam;

uitableTracers = uitable('parent', hFigAddLagTracer,...
    'data', TracerInit,...
    'position', [10 90 400 30]);

myTracerhandles=guidata(hFigAddLagTracer);

WidthMenus = uibuttonparam.Position(1)+ uibuttonparam.Extent(3);
HeightInit = uitableTracers.Position(2)+uitableTracers.Position(4);
hFigAddLagTracer.Position(3) = WidthMenus + pad;
hFigAddLagTracer.Position(4) = HeightInit + pad;

hApplyTracer.Position(3) = hApplyTracer.Extent(3) + pad;
hApplyTracer.Position(4) = hApplyTracer.Extent(4) + pad;
hApplyTracer.Position(1) = hFigAddLagTracer.Position(3) - (hApplyTracer.Position(3) + pad);

hApplyTracer.Callback = @ApplyTracersCallback;

uibuttoncoord.Callback = @GetTracerCoord;

function GetTracerCoord(~,~)
    myTracerhandles.TracerNumS = uiTracerNum.String;
    myTracerhandles.FrameNumS = uiFrameNum.String;    
end

function ApplyTracersCallback(~,~)
        save('TracerInit.mat', 'TracerInit')
end

function AddParam(~,~)
    
    myTracerhandles.FrameNumS = uiFrameNum.String;
    myTracerhandles.PopMenS = uipopmenuparam.String;
    myTracerhandles.PopMenV = uipopmenuparam.Value;
    myTracerhandles.TracerNumS = uiTracerNum.String;
    myTracerhandles.TracerInit = TracerInit;
    
    VAR_Tracer = str2double(myTracerhandles.TracerNumS);
    VAR_Frame = str2double(myTracerhandles.FrameNumS);
    VAR_Param_Name = char(myTracerhandles.PopMenS(myTracerhandles.PopMenV));
    
    TracerInit{VAR_Tracer, 1} = VAR_Tracer;
    TracerInit{VAR_Tracer, 2} = VAR_Frame;
    Line = TracerInit(VAR_Tracer,:);
    S = size(Line);
    
    POS = 5; % first parameter in position 5
    if S(2) < POS
        % no parameter entered yet
        TracerInit{VAR_Tracer, POS} = VAR_Param_Name;
    else
        TF = isempty(TracerInit{VAR_Tracer, POS});
        
        while(TF == 0) 
            POS = POS+1;
            if S(2) < POS
                TF = 1;
            else
                TF = isempty(TracerInit{VAR_Tracer, POS});
            end
        end
        TracerInit{VAR_Tracer, POS} = VAR_Param_Name;
    end
    
    uitableTracers.Data = TracerInit;
    uitableTracers.Position(3) =  uitableTracers.Extent(3);
    uitableTracers.Position(4) =  uitableTracers.Extent(4);
    
    WidthTable = uitableTracers.Position(1)+uitableTracers.Position(3);
    WidthMenus = uibuttonparam.Position(1)+ uibuttonparam.Extent(3);
    
    if WidthMenus > WidthTable
        hFigAddLagTracer.Position(3) = WidthMenus + pad;
    else
        hFigAddLagTracer.Position(3) = WidthTable + pad;
    end
    
    HeightInit = uitableTracers.Position(2)+uitableTracers.Position(4);
    hFigAddLagTracer.Position(4) = HeightInit + pad;
    hApplyTracer.Position(1) = hFigAddLagTracer.Position(3) - (hApplyTracer.Position(3) + pad);
         
end


end

