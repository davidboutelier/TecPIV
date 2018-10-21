% script seting up vector stats


%% inputs to be passed to the function
F=99; % frame number
Plottype = 'density'; %'UVscatter'; % type of plot 'histogram', 'subpixel'
exppath=fullfile('/Users','david','Documents','PIV_STAGE','SS2');
vectorpath=fullfile('Raw','Vectors');

load(fullfile(exppath,vectorpath,['Vector_',num2str(F),'.mat']));

S=size(U);
figure(1)

switch Plottype
    case 'histogram'
        h1=histogram(U)
        h1.BinWidth = 0.1;
        hold on
        h2=histogram(V)
        h2.BinWidth = 0.1;
        legend('U', 'V','location','northwest');
        xlabel('pixel');
    
    case 'subpixel'
        h3=histogram(U-fix(U))
        h3.BinWidth = 0.1;
        hold on
        h4=histogram(V-fix(V))
        h4.BinWidth = 0.1;
        legend('U-fix(U)', 'V-fix(V)','location','northwest')
        xlabel('pixel');
    
    case 'UVscatter'
        plot(reshape(U,[],1),reshape(V,[],1),'+')
        xlabel('U (pixel)');
        ylabel('V (pixel)');
        
    case 'density'
        [values, centers] = hist3([reshape(U,[],1) reshape(V,[],1)],[101 101]);
        imagesc(centers{:}, values.')
        colorbar
        axis xy
        
        
end
        








