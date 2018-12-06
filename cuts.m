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
