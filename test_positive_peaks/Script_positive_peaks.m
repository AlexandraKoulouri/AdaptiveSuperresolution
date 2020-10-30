% This is the test case appeared in Adaptive Superresolution in
% Deconvolution of sparse peaks for the cases of 5 positive peaks 

close all; clear
addpath(genpath(cd)) 

% Define Peak locations and Amplitudes 
xi = [0.49 0.486 0.4 0.613 0.41];
yi =  [0.56 0.65 0.47 0.5  0.44];

c  = [1.2 0.8 1 1 1.5];

%Define observation grid
m_size = 40;
um = linspace(0,1,m_size);
[x,y] = meshgrid(um,um);
xm =  x(:);
ym =  y(:);
clear x y
sigma = [2*um(2) ;2.5*um(2)]; %always the kernel size wider than the pixel size! %[2*um(2);2.5*um(2)];
alpha = [1; 0];
ob= SimulateObservationsGaussian(xm,ym,sigma,c,xi,yi,alpha);
n = randn(length(ob),1);
 
SNR = 40; %Noise level
scale = (norm(ob)^2/norm(n)^2)* 10^(-SNR/10); 
ob_n = ob' + sqrt(scale)*n;
ob_n(ob_n<1e-10) = 0;

%% Adaptive algorithm: Define the initial computational grid of the adaptive method
Lim = 0.01*max(sigma);
size_comp_grid = 15;
QuietPlot = 0;
[xi_es,c_es]= AS_alg(ob_n,xm,ym,sigma,alpha,size_comp_grid,Lim,QuietPlot,xi,yi);

%% Plot final result!
figure;
set(gcf, 'Units','centimeters', 'Position',[1 1 14 7]) 

%PLOT LOW_RESOLUTION IMAGE
axes('position',[0.01 0.2 0.3 0.65])
u = linspace(0,1,100);
[xc yc] = meshgrid(u,u);
xr = xc(:);
yr = yc(:);
M_low = InterpolateMatrix2d(delaunay(xm,ym),[xm ym],[xr(:) yr(:)]);
ob_up = full(M_low*ob_n);
LowResImg = reshape(ob_up,length(u),length(u));
imagesc(LowResImg)
title('Low resolution image','fontsize',9)
colormap('hot')
axis square
axis off
axis xy
set(gca,'fontsize',7)
h_axis = axis;
%axis(axis_h./m_size)
set(gca,'Xtick',[h_axis(1)+0.5,h_axis(2)],'XTickLabel',{'0','1'})
set(gca,'Ytick',[h_axis(3),h_axis(4)],'YTickLabel',{'0','1'})



% THIS IS THE SUPERRESOLUTION RESULT
axes('position',[0.34 0.2 0.3 0.65])
u = linspace(0,1,340);
[x y] = meshgrid(u,u);
xr = x(:);
yr = y(:);
M_rec = InterpolateMatrix2d(delaunay(xr,yr),[xr yr],[xi_es(:,1) xi_es(:,2)]);
c_es_grid = full(M_rec'*c_es);
SuperResImg = reshape(c_es_grid,length(u),length(u));
h = fspecial('gaussian', 10,1.3);
SuperResImg = imfilter(SuperResImg,h);
imagesc(SuperResImg)
axis square
axis xy
axis off
set(gca,'fontsize',7)
set(gca,'Xtick',[h_axis(1),h_axis(2)],'XTickLabel',{'0','1'})
set(gca,'Ytick',[h_axis(3),h_axis(4)],'YTickLabel',{'0','1'})
title('Super-resolution Result','fontsize',8)

axes('position',[0.35 0.12 0.3 0.05])
colorbar('location','South')
colormap('hot')
caxis([min(c_es) max(c_es)])  
axis off


% Plot the locations of the actual and reconstucted peaks
axes('position',[0.685 0.2 0.3 0.65])

for ii=1:length(c_es)
    if c_es(ii)>0
          plot(xi_es(ii,1),xi_es(ii,2),'o','markersize',6,'MarkerFaceColor',[0.7 0.7 0.7])%,'MarkerEdgeColor','none')
         hold on
    else
         plot(xi_es(ii,1),xi_es(ii,2),'o','markersize',6,'linewidth',1)%,'MarkerFaceColor',[55  150 171]./255)%,'MarkerEdgeColor','none')
         hold on
    end
end
for ii = 1:length(c)
    if c(ii)>0
       plot(xi(ii),yi(ii),'x','markersize',10,'color',[0.7 0.1 0.1])%,'MarkerEdgeColor','none')
    else
       plot(xi(ii),yi(ii),'+','markersize',10,'color',[55  150 171]./255,'linewidth',1)%,'MarkerEdgeColor','none')            `
    end
end
axis([min(xi_es(:,1))-0.1*min(xi_es(:,1))   max(xi_es(:,1))+0.1*min(xi_es(:,1))    min(xi_es(:,2))-0.1*min(xi_es(:,2))    max(xi_es(:,2))+0.1*min(xi_es(:,2))])

title('Peak Locations','FontSize',8)
axis square
axis xy
box on
axis on

%axis off
%axes('position',[0.05 0.15 0.40 0.65])
axes('position',[0.80 0.041 0.12 0.05])
plot(3,0,'x','color',[0.7 0.091 0.1]);
hold on
plot(3,-0.2,'o','markersize',6,'MarkerFaceColor',[0.7 0.1 0.1])
hold on

text(-1.1,0.0,'Original peaks','FontSize',10)
text(-1.1,-0.2,'Estimated peaks','FontSize',10)

axis off



%%Similary measure!
[SLE,Flow] = SigSimularityEMD(c_es,c(:),xi_es,[xi yi]);
