function DisplayIMG(Ax,F, Datasets, Value,RawCpt,VectorField,Derivative)
%This function makes the figures. The image is not passed on to the function, only the dataset and frame number 
% F: the frame number to be displayed
% Datasets: the cell array with all the ddatasets in the current project
% Value: the number of the dataset to be employed
% RawCpt: a cell array with the colormaps
% VectorField: a cell array with the attributs to be used for the vectors
% Derivative: a cell array with the attributs to be used for the derivatives


    % Commun parameters regardless of data type
    TecPivFolder =  RawCpt{1,7};
    PIVStagePath=Datasets{Value,2}; % where the project was created (in the computer)
    ProjectName=Datasets{Value,3}; % name of project
    Dt = Datasets{Value,7}; % time increment for the dataset
    ValueDataImg = Datasets{Value,11}; % number of the dataset containing the images
    ImgPath = Datasets{ValueDataImg,1}; % path of the images
    DatasetName = Datasets{Value,1};

    % get the datatype from the datasets
    DataType=Datasets{Value,15};
    
        
    switch DataType
        case 'image'
            % Plot image
            imgFramepath=fullfile(PIVStagePath,ProjectName,ImgPath,['IMG_' num2str(F) '.tif']);
            I=imread(imgFramepath);
    
            TecPIV_Display_Image(I,Ax,RawCpt,VectorField,Derivative)
            
        case 'vector'
            % get background image of model
            imgFramepath=fullfile(PIVStagePath,ProjectName,ImgPath,['IMG_' num2str(F) '.tif']);
            I=imread(imgFramepath);
            
            MaskDatasetNumber = Datasets{Value,17}; % dataset number where mask is
            MaskPath = Datasets{MaskDatasetNumber,1}; % path of the mask
                
            DATA = load(fullfile(PIVStagePath,ProjectName,MaskPath,['Mask_',num2str(F),'.mat']));
            ROI = DATA.ROI;
            RoiMask = DATA.RoiMask;
            clear DATA
            
            MaskExist = Datasets{Value,16};
            
%             % check if there is a mask
%             MaskExist = Datasets{Value,16};
%             if MaskExist == 1 % there is a mask
%                 MaskDatasetNumber = Datasets{Value,17}; % dataset number where mask is
%                 MaskPath = Datasets{MaskDatasetNumber,1}; % path of the mask
%                 
%                 DATA = load(fullfile(PIVStagePath,ProjectName,MaskPath,['Mask_',num2str(F),'.mat']));
%                 ROI = DATA.ROI;
%                 RoiMask = DATA.RoiMask;
%                 clear DATA
%                 
%             else
%                 ROI = [];
%                 RoiMask = [];
%             end 
            
            RangeType=Derivative{1,3};
            DeriveCpt=Derivative{1,4};
            
             ListMATLABCPT={'parula','jet','hsv','hot','cool','spring','summer',...
                'autumn','winter','gray','bone','copper','pink','lines',...
                'colorcube','prism','flag','white'};
             Lia = ismember(DeriveCpt,ListMATLABCPT); % check if colormap is in the list of MATLAB colormaps
                
             if Lia == 1 % it is a MATLAB colormap
                    RGB=DeriveCpt;
                    %colormap(gca,RGB);
             else % it is a cutom colormap loaded from the folder
                    load(fullfile(TecPivFolder,'toolbox','colormaps',DeriveCpt));
                    %colormap(gca,RGB);
             end 
                         
             TecPIV_Display_Vectors(I,Ax,RawCpt,VectorField,Derivative,Dt,MaskExist,ROI,RoiMask,RGB,RangeType)
        
        case 'Lagrangian'
            
            % Plot image
            imgFramepath=fullfile(PIVStagePath,ProjectName,ImgPath,['IMG_' num2str(F) '.tif']);
            I=imread(imgFramepath);
            
            % check if there is a mask
            MaskExist = Datasets{Value,16};
            if MaskExist == 1 % there is a mask
                MaskDatasetNumber = Datasets{Value,17}; % dataset number where mask is
                MaskPath = Datasets{MaskDatasetNumber,1}; % path of the mask
                DATA = load(fullfile(PIVStagePath,ProjectName,MaskPath,['Mask_',num2str(F),'.mat']));
                ROI = DATA.ROI;
                RoiMask = DATA.RoiMask;
                clear DATA
                
            else
                ROI = [];
                RoiMask = [];
            end
            
            % laod the lagrangian sum vectors
            DATA = load(fullfile(PIVStagePath,ProjectName, DatasetName, ['Vector_',num2str(F),'.mat']));
            X = DATA.X;
            Y = DATA.Y;
            U = DATA.U;
            V = DATA.V;
            clear DATA
            
            DATA = load(fullfile(PIVStagePath,ProjectName,'Datasets.mat'));
            SaveDatasets = DATA.SaveDatasets;
            clear DATA
            
            SourceSumNumberDataset = SaveDatasets{Value,18};
            SourceSumDataName = SaveDatasets{SourceSumNumberDataset,1};
            
            DATA = load(fullfile(PIVStagePath,ProjectName, 'PIV_parameters.mat'));
            SavePIVparam = DATA.SavePIVparam;
            
            NumPass = SavePIVparam{11,1};
            
            switch (NumPass)
                case 1
                    DX = SavePIVparam{1,1};
                case 2
                    DX = SavePIVparam{12,1};
                case 3
                    DX = SavePIVparam{13,1};
                case 4
                    DX = SavePIVparam{14,1};
            end
            % check if there is a mask
            MaskExist = Datasets{Value,16};
            if MaskExist == 1 % there is a mask
                MaskDatasetNumber = Datasets{Value,17}; % dataset number where mask is
                MaskPath = Datasets{MaskDatasetNumber,1}; % path of the mask
                
                DATA = load(fullfile(PIVStagePath,ProjectName,MaskPath,['Mask_',num2str(F),'.mat']));
                ROI = DATA.ROI;
                RoiMask = DATA.RoiMask;
                clear DATA
                
            else
                ROI = [];
                RoiMask = [];
            end 
            
            RangeType=Derivative{1,3};
            DeriveCpt=Derivative{1,4};
            
             ListMATLABCPT={'parula','jet','hsv','hot','cool','spring','summer',...
                'autumn','winter','gray','bone','copper','pink','lines',...
                'colorcube','prism','flag','white'};
             Lia = ismember(DeriveCpt,ListMATLABCPT); % check if colormap is in the list of MATLAB colormaps
                
             if Lia == 1 % it is a MATLAB colormap
                    RGB=DeriveCpt;
                    %colormap(gca,RGB);
             else % it is a cutom colormap loaded from the folder
                    load(fullfile(TecPivFolder,'toolbox','colormaps',DeriveCpt));
                    %colormap(gca,RGB);
             end 
             
            TecPIV_Display_LagGrid(I,Ax,RawCpt,VectorField,Derivative,Dt,MaskExist,ROI,RoiMask,RGB,RangeType, X,Y,U,V,DX)
            
    
    end
end

          
            
%         case 'Lagrangian'
%             % load images
%             imgFramepath=fullfile(PIVStagePath,ProjectName,ImgPath,['IMG_' num2str(F) '.tif']);
%             I=imread(imgFramepath);
%     
%             % adjust image to requested range
%             I = imadjust(I,[RelativeMin RelativeMax],[0 1]);
%             
%             % laod the lagrangian sum vectors
%             load(fullfile(PIVStagePath,ProjectName, DatasetName, ['Vector_',num2str(F),'.mat']));
%             
%             %load the masks
%             MaskDatasetNumber = Datasets{Value,17}; % dataset number where mask is
%             MaskPath = Datasets{MaskDatasetNumber,1};
%             load(fullfile(PIVStagePath,ProjectName,MaskPath,['Mask_',num2str(F),'.mat']));
%             
%             % patch the holes in grids
%             XX=inpaint_nans(X,3);
%             YY=inpaint_nans(Y,3);
%             UU=inpaint_nans(U,3);
%             VV=inpaint_nans(V,3);
%             
%             % get the minmax of vector field
%             c = 0;
%             if c == 1
%                 xmin=0;
%                 xmax=size(I,2)-1;
%                 ymin=0;
%                 ymax=size(I,1)-1;
%             else
%                 xmin = floor(min2(XX));
%                 ymin = floor(min2(YY));
%                 xmax = ceil(max2(XX));
%                 ymax = ceil(max2(YY));
%             end
%             
%             xwidth=[xmin xmax];
%             ywidth=[ymin ymax];
%             
%             DX=abs(XX(1,1)-XX(1,2));
%             DY=abs(YY(1,1)-YY(2,1));
%             
%             % define new x and y vectors
%             Nx=xmin:1:xmax;
%             Ny=ymin:1:ymax;
% 
%             % create new x and y grids
%             [NX,NY]=meshgrid(Nx,Ny);
%             
%             % calculate gradient of lagrangian displacements
%             % points are initially every 64 pixels (PIV resolution for each experiment)
%             [exx,exy] = gradient(UU,DX,DY);
%             [eyx,eyy] = gradient(VV,DX,DY);
%             
%             % interpolate the gradients to get finer resolution
%             Nexx = griddata(XX,YY,exx,NX,NY,'linear');
%             Nexy = griddata(XX,YY,exy,NX,NY,'linear');
%             Neyy = griddata(XX,YY,eyy,NX,NY,'linear');
%             Neyx = griddata(XX,YY,eyx,NX,NY,'linear');
%             
%             Derive = 0.5*(exy+eyx +(exx.*exy+eyx.*eyy));
%             NDerive = 0.5*(Nexy+Neyx+(Nexx.*Nexy+Neyx.*Neyy));
%             
%             
%             % check that mask exist
%             MaskExist = Datasets{Value,16};
%             
%             if MaskExist == 1 % if mask does exist
%                 load(fullfile(PIVStagePath,ProjectName,MaskPath,['Mask_',num2str(F),'.mat']));
%                 XM=RoiMask(:,1);
%                 YM=RoiMask(:,2);
%                 mask = poly2mask(XM-xmin,YM-ymin,size(NDerive,1),size(NDerive,2));
%                 
%                 % delete extrapolated outside the mask area
%                 NDerive(mask == 0) = 0;
%             end
%             
% 
%             MaxRange=max(max(Derive));
%             MinRange=min(min(Derive));
%             Range=[MinRange, MaxRange];
%             
%             DerivROI=imref2d(size(NDerive),xwidth,ywidth);
%             
%             nx=4; ny=4;
%             XS = X(1:nx:end, 1:ny:end);
%             YS = Y(1:nx:end, 1:ny:end);
%             US = U(1:nx:end, 1:ny:end);
%             VS = V(1:nx:end, 1:ny:end);
%             %typevectorS = typevector(1:nx:end, 1:ny:end);
% 
%             XS=XS(1,:);
%             YS=(YS(:,1))';
%             
%             I2=im2double(I); %convert image from uitnt16 to double
%             IRGB = double2rgb(I2, gray); % convert to RGB
%             Ider=mat2im(NDerive,parula,Range); % convert matrix to RGB
%             
%             RA=imref2d(size(I)); % area of image one (background)
%             [D,RD] = imfuse(IRGB,RA,Ider,DerivROI,'method','blend'); % blend the two images
%             
%             subimage(D);
%             iDeriv2=colorbar('location','westoutside');
%             set(get(iDeriv2,'child'),'YData',Range);
%             set(iDeriv2,'YLim',Range);  
%             caxis(Range);
%             set(gca, 'TickDir', 'out')
%             
%             % calculate the outline of the interpolated grid
%             k=1;
%             XX2=XX(1:k:end,1:k:end);
%             YY2=YY(1:k:end,1:k:end);
%             
%             N = size(XX2,1);
%             M = size(XX2,2);
%             
%             for i=1:M-1
%                 for j=1:N
%                     VX = [XX2(j,i), XX2(j,i+1)];
%                     VY = [YY2(j,i), YY2(j,i+1)];
%                     plot(VX,VY,'-','color','black')
%                     hold on
%                 end
%             end
% 
%             % then draw vertical lines
%             for i=1:M
%                 for j=1:N-1
%                     VX = [XX2(j,i), XX2(j+1,i)];
%                     VY = [YY2(j,i), YY2(j+1,i)];
%                     plot(VX,VY,'-','color','black')
%                     hold on
%                 end
%             end
%             daspect([1 1 1])
%             
%     
%     end
%         
          

    
    


    
    
            

            


    
    
                 


