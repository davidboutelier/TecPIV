clc
close all
clear all
%Type='Divergent-Linear-Lab';
Type='Divergent-Parametric-Lab';

StartColor='Dark_Blue';
CenterColor='White';
EndColor='Red';

Ncolors=256;

L=[0:1:100];
a=[-100:1:100];
b=a;


switch StartColor
    case 'Blue'
        Lab=rgb2lab([0 0.4 0.8]);
    case 'Dark_Blue'
        Lab=rgb2lab([0 0.2 0.4]);
    case 'Light_Blue'
        Lab=rgb2lab([0.8 0.7 1]); 
    case 'White'
        Lab=rgb2lab([1 1 1]);
    case 'Dark_Red'
        Lab=rgb2lab([0.4 0 0]);
    case 'Red'
        Lab=rgb2lab([0.8 0 0]);
    case 'Dark_Grey'
        Lab=rgb2lab([0.25 0.25 0.25]);
    case 'Light_Grey'
        Lab=rgb2lab([0.9 0.9 0.9]);
    case 'Grey'
        Lab=rgb2lab([0.5 0.5 0.5]);
       
end
SL=Lab(1);
Sa=Lab(2);
Sb=Lab(3);

switch EndColor
    case 'Blue'
        Lab=rgb2lab([0 0.4 0.8]);

    case 'Dark_Blue'
        Lab=rgb2lab([0 0.2 0.4]);

    case 'Light_Blue'
        Lab=rgb2lab([0.8 0.7 1]);
        
    case 'White'
        Lab=rgb2lab([1 1 1]);
    case 'Dark_Red'
        Lab=rgb2lab([0.4 0 0]);
    case 'Red'
        Lab=rgb2lab([0.8 0 0]);
        
    case 'Dark_Grey'
        Lab=rgb2lab([0.25 0.25 0.25]);
    case 'Light_Grey'
        Lab=rgb2lab([0.9 0.9 0.9]);
    case 'Grey'
        Lab=rgb2lab([0.5 0.5 0.5]);

end
EL=Lab(1);
Ea=Lab(2);
Eb=Lab(3);

switch CenterColor
    case 'White'
        Lab=rgb2lab([1 1 1]);
        case 'Dark_Grey'
        Lab=rgb2lab([0.25 0.25 0.25]);
    case 'Light_Grey'
        Lab=rgb2lab([0.9 0.9 0.9]);
    case 'Grey'
        Lab=rgb2lab([0.5 0.5 0.5]);
end

CL=Lab(1);
Ca=Lab(2);
Cb=Lab(3);

switch Type
    case 'Linear-Lab'
        x=EL-SL;
        y=Ea-Sa;
        z=Eb-Sb;
        [azimuth,elevation,r] = cart2sph(x,y,z);
        
        for i=1:Ncolors-2
            [X(i),Y(i),Z(i)] = sph2cart(azimuth,elevation,i*(r/(Ncolors-1)));
            X(i)=X(i)+SL;
            Y(i)=Y(i)+Sa;
            Z(i)=Z(i)+Sb;
        end
        
        X=[SL,X,EL];
        Y=[Sa,Y,Ea];
        Z=[Sb,Z,Eb];
    
    case 'Divergent-Linear-Lab'
        
        if mod(Ncolors,2) == 0
            disp('W: number of colors changed to nearest odd number')
            Ncolors=Ncolors+1;
        end
        x=CL-SL;
        y=Ca-Sa;
        z=Cb-Sb;
        [azimuth,elevation,r] = cart2sph(x,y,z);
        
        for i=1:(Ncolors-3)/2
            [Xa(i),Ya(i),Za(i)] = sph2cart(azimuth,elevation,i*(r/((Ncolors-3)/2 +1)));
            Xa(i)=Xa(i)+SL;
            Ya(i)=Ya(i)+Sa;
            Za(i)=Za(i)+Sb;
        end
        
        x=EL-CL;
        y=Ea-Ca;
        z=Eb-Cb;
        [azimuth,elevation,r] = cart2sph(x,y,z);
        
        for i=1:(Ncolors-3)/2
            [Xb(i),Yb(i),Zb(i)] = sph2cart(azimuth,elevation,i*(r/((Ncolors-3)/2 +1)));
            Xb(i)=Xb(i)+CL;
            Yb(i)=Yb(i)+Ca;
            Zb(i)=Zb(i)+Cb;
        end
        
        X=[SL,Xa,CL,Xb,EL];
        Y=[Sa,Ya,Ca,Yb,Ea];
        Z=[Sb,Za,Cb,Zb,Eb];
        
        
        
    
    case 'Divergent-Parametric-Lab'
        Xi=[SL,CL,EL];
        Yi=[Sa,Ca,Ea];
        Zi=[Sb,Cb,Eb];
        pt = interparc(Ncolors,Xi,Yi,Zi);
        
        X=pt(:,1)';
        Y=pt(:,2)';
        Z=pt(:,3)';

        
        
end



RGB=lab2rgb([X' Y' Z']);

RGB(RGB<0)=0;
RGB(RGB>1)=1;






figure(1)
s=24;
scatter3(X,Y,Z,s)
xlabel('L*')
ylabel('a*')
zlabel('b*')
colormap(RGB)
colorbar

figure(2)
colormap(RGB)
imagesc(sineramp);
axis equal tight
colorbar


