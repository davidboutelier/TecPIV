function TecPIV_Lag_Sum(DataSets,ThisDataSetNumber)
% This function calculate the lagrangian sum of incremental displacements
%
%% Parameters from GUI
DatasetFolder = DataSets{ThisDataSetNumber,1};
PathData = DataSets{ThisDataSetNumber,2};
ProjectID = DataSets{ThisDataSetNumber,3};
        
NumberImages = DataSets{ThisDataSetNumber,4};
ImageWidth = DataSets{2,5};
ImageHeight = DataSets{2,6};
TimeInc = DataSets{ThisDataSetNumber,7};
ImageInc = DataSets{ThisDataSetNumber,8};
StartNumber = DataSets{ThisDataSetNumber,9}; 
EndNumber = DataSets{ThisDataSetNumber,10}; 

ImageFolderNumber = DataSets{ThisDataSetNumber,11}; % number of the dataset that contains the background images;

%% Parameters for summation
Npt=9;
Nsteps=(EndNumber-StartNumber)/ImageInc + 1;

%% create folder if it doesn't exist already.
Framepath=fullfile(PathData,ProjectID,DatasetFolder);
cd(Framepath);

NumberDatasets = size(DataSets,1);
NLag = 0;
for j = 2:NumberDatasets
    DatasetName =  DataSets{j,1}
    TF = contains(DatasetName,'Lagrangian');
    if TF == 1
        NLag = NLag+1;
    end  
end

NewFolderName=['Lagrangian_Sum_' num2str(NLag+1)]; 
mkdir(NewFolderName);


% 
% if ~exist('Lagrangian_Sum', 'dir')
%   mkdir('Lagrangian_Sum');
% end
cd(fullfile(PathData,ProjectID));

% %% Display progress
% progressStepSize = 1;
% ppm = ParforProgMon('Summing vector fields:  ', Nsteps, progressStepSize, 500, 40);

obj = ProgressBar(Nsteps,'Title','Summing vectors');

%% Main loop
for j=1:Nsteps
    F=StartNumber+(j-1)*ImageInc; % Framenumber
    
    % load vector file and get size of vector field
    load(fullfile(Framepath, ['Vector_' num2str(F) '.mat']),'X','Y','U','V');
    [Nx,Ny]=size(X);
    N(j)=Nx*Ny;
    
    if j == 1
        Nxi=Nx;
        Nyi=Ny;
        
        CX=X;  
        CY=Y;
        CU=U;
        CV=V;
        
        D = TimeInc.*ImageInc.*(U.^2+V.^2).^0.5;
        CD = D;
        
        % define target position (X,Y)
        Xt=X+U; 
        Yt=Y+V;

        Xt=reshape(Xt,N(1),1); L = length(Xt);
        Yt=reshape(Yt,N(1),1);
        
        clear X Y U V D
        X=CX;
        Y=CY;
        U=CU;
        V=CV;
        D=CD;
        
        SaveName=fullfile(Framepath, NewFolderName,['Vector_' num2str(F) '.mat']);
        save(SaveName,'X','Y','U','V','D'); 
        clear X Y U V D
    else
        X=reshape(X,N(j),1);
        Y=reshape(Y,N(j),1);
        
        U=reshape(U,N(j),1);
        V=reshape(V,N(j),1);
        
        % preallocate the target disp
        Ut=zeros(length(Xt),1);
        Vt=zeros(length(Xt),1);
        %indSort = zeros(Npt,1);
        
        % filter with image size
        Ut(Xt >= ImageWidth) = 0;
        Ut(Xt <= 0) = 0;
        Vt(Yt >= ImageHeight) = 0;
        Vt(Yt <= 0) = 0;
        
        Xt(Xt >= ImageWidth) = ImageWidth;
        Xt(Xt<= 0) = 0;
        Yt(Yt>= ImageHeight) = ImageHeight;
        Yt(Yt<= 0) = 0;
        
        parfor i=1:L
            TX=Xt(i);
            TY=Yt(i);

            if ~isnan(TX) && ~isnan(TY)
                
                if (TX >= ImageWidth) || (TX <= 0) || (TY >= ImageHeight) || (TY <= 0)
                    if TX >= ImageWidth
                        Ut(i) = 0;
                        Xt(i) = ImageWidth;
                    end
                    if TX <= 0 
                        Ut(i) = 0;
                        Xt(i) = 0;
                    end
                    if TY >= ImageHeight
                        Vt(i) = 0;
                        Yt(i) = ImageHeight;
                    end
                    if TY <= 0
                        Vt(i) =0;
                        Yt(i) = 0;
                    end
                else

                    % compute Euclidean distances:
                    dist=bsxfun(@hypot,Y-TY,X-TX);

                    % find Npt minimum distances
                    indSort=bsxfun(@TecPIV_NMin,dist,Npt);

                    Xsub=X(indSort);
                    Ysub=Y(indSort);
                    Usub=U(indSort);
                    Vsub=V(indSort);

                    % check if first 3 are colinear
                    k=1;
                    tf =1;

                    while (tf~=0 && k<=Npt-2)
                        p1 = [Xsub(1) Ysub(1)];
                        p2 = [Xsub(2) Ysub(2)];
                        p3 = [Xsub(k+2) Ysub(k+2)];
                        tf = TecPIV_collinear(p1,p2,p3);
                        k=k+1;
                    end
                
                    if k == Npt-1
                        Ut(i)=TecPIV_IDW(Xsub,Ysub,Usub,TX,TY,-2,'ng',Npt);
                        Vt(i)=TecPIV_IDW(Xsub,Ysub,Vsub,TX,TY,-2,'ng',Npt);
                        disp('error, could not find 3 non-colinear points, fallback to inverse distance weighting method')
                    else
                        % calculate the equation of the plane
                        A = [Xsub(1) Ysub(1) Usub(1)];
                        B = [Xsub(2) Ysub(2) Usub(2)];
                        C = [Xsub(k+1) Ysub(k+1) Usub(k+1)];

                        AB = [B(1)-A(1) B(2)-A(2) B(3)-A(3)];
                        AC = [C(1)-A(1) C(2)-A(2) C(3)-A(3)];

                        s1 = AB(2)*AC(3)-AB(3)*AC(2);
                        s2 = AB(3)*AC(1)-AB(1)*AC(3);
                        s3 = AB(1)*AC(2)-AB(2)*AC(1);

                        s4 = -1*(s1*A(1)+s2*A(2)+s3*A(3));
                        UT = (-s4 -(s1*TX+s2*TY))/s3;

                        A = [Xsub(1) Ysub(1) Vsub(1)];
                        B = [Xsub(2) Ysub(2) Vsub(2)];
                        C = [Xsub(k+1) Ysub(k+1) Vsub(k+1)];

                        AB = [B(1)-A(1) B(2)-A(2) B(3)-A(3)];
                        AC = [C(1)-A(1) C(2)-A(2) C(3)-A(3)];

                        s1 = AB(2)*AC(3)-AB(3)*AC(2);
                        s2 = AB(3)*AC(1)-AB(1)*AC(3);
                        s3 = AB(1)*AC(2)-AB(2)*AC(1);

                        s4 = -1*(s1*A(1)+s2*A(2)+s3*A(3));
                        VT = (-s4 -(s1*TX+s2*TY))/s3;

                        Ut(i) = UT;
                        Vt(i) = VT;
                    end
                end
            end
        end
        
        % write output
        
        U = reshape(Ut,Nxi,Nyi);
        V = reshape(Vt,Nxi,Nyi);
        D = TimeInc.*ImageInc.*(U.^2+V.^2).^0.5;
        
        CU = CU + U;
        CV = CV + V;
        CX = reshape(Xt,Nxi,Nyi);
        CY = reshape(Yt,Nxi,Nyi);
        
        CD = CD + D;
        
        clear X Y U V D
        X=CX;
        Y=CY;
        U=CU;
        V=CV;
        D=CD;
        
        Xt = Xt + Ut;
        Yt = Yt + Vt;
        
        SaveName=fullfile(Framepath, NewFolderName,['Vector_' num2str(F) '.mat']);
        save(SaveName,'X','Y','U','V', 'D');
  
    end
  

obj.step([], [], []); 
end

obj.release();




