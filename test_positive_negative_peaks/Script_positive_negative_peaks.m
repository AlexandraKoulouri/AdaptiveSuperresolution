close all; clear
addpath(genpath(cd)) 
% the signals are positive and negative
close all; clear

% Define Peak locations and Amplitudes 
xi = [0.49 0.486 0.4 0.613];
yi =  [0.56 0.65 0.47 0.5];

c  = [-1 0.8 -1 1];

%Define observation grid and simulate observation
m_size = 40;
um = linspace(0,1,m_size);
[x,y] = meshgrid(um,um);
xm =  x(:);
ym =  y(:);
clear x y
sigma = [2*um(2) ;2.5*um(2)]; %always the kernel size wider than the pixel size! %[2*um(2);2.5*um(2)];
alpha = [0.8; 0.2];
ob= SimulateObservationsGaussian(xm,ym,sigma,c,xi,yi,alpha);
n = randn(length(ob),1);

SNR =40;
scale = (norm(ob)^2/norm(n)^2)* 10^(-SNR/10); 
ob_n = ob' + sqrt(scale)*n;


%% Adaptive algorithm: Define the initial computational grid of the adaptive method
Lim = 0.2*um(2);%max(sigma);
size_comp_grid = 15;
QuietPlot = 1;
[xi_es,c_es]= AS_alg(ob_n,xm,ym,sigma,alpha,size_comp_grid,Lim,QuietPlot,[],[]);

%% Plot final result!
figure
set(gcf, 'Units','centimeters', 'Position',[1 1 16 7]) 
axes('position',[0.045 0.22 0.29 0.65])
%observation image
%u = linspace(0,1,140);
%[x y] = meshgrid(u,u);
%xr = x(:);
%yr = y(:);
%M_ob = InterpolateMatrix2d(delaunay(xm,ym),[xm ym],[xr yr]);
%ob_1 = M_ob*ob';
%obImg = reshape(ob_1,length(u),length(u));
obImg = reshape(ob,m_size,m_size);
imagesc(obImg)
%title('Observations','fontsize',8)
%h_axis = axis;
axis([0.3 0.7 0.35 0.75]*size(obImg,1))
h_axis= [0.3 0.7 0.35 0.75]*size(obImg,1);
set(gca,'Xtick',[h_axis(1),h_axis(2)],'XTickLabel',{'0.3','0.7'})
set(gca,'Ytick',[h_axis(3),h_axis(4)],'YTickLabel',{'0.35','0.75'})
axis square
axis xy
set(gca,'fontsize',9)
hold on
xi_pix = 1+xi*(size(obImg,1)-1);
yi_pix = 1+yi*(size(obImg,2)-1);
for ii = 1:length(c)
    if c(ii)>0
       plot(xi_pix(ii),yi_pix(ii),'x','markersize',10,'color',[0.7 0.1 0.1])%,'MarkerEdgeColor','none')
    else
       plot(xi_pix(ii),yi_pix(ii),'+','markersize',10,'color',[55  150 171]./255,'linewidth',1)%,'MarkerEdgeColor','none')            `
    end
end

box on
axes('position',[0.375 0.22 0.3 0.65])
u = linspace(0,1,120);
[x y] = meshgrid(u,u);
xr = x(:);
yr = y(:);
M_rec = InterpolateMatrix2d(delaunay(xr,yr),[xr yr],[xi_es(:,1) xi_es(:,2)]);
c_es_grid = full(M_rec'*c_es);
SuperResImg = reshape(c_es_grid,length(u),length(u));
%h = fspecial('gaussian', 5,5);
h = fspecial('disk',2)
SuperResImg = imfilter(SuperResImg,h);
imagesc(SuperResImg)
set(gca,'fontsize',8)
axis square
axis xy
hold on
xi_pix = 1+xi*(size(SuperResImg,1)-1);
yi_pix = 1+yi*(size(SuperResImg,2)-1);
for ii = 1:length(c)
    if c(ii)>0
       plot(xi_pix(ii),yi_pix(ii),'x','markersize',10,'color',[0.7 0.1 0.1])%,'MarkerEdgeColor','none')
    else
       plot(xi_pix(ii),yi_pix(ii),'+','markersize',10,'color',[55  150 171]./255,'linewidth',1)%,'MarkerEdgeColor','none')            `
    end
end



%set(gca,'fontsize',7)
axis([0.3 0.7 0.35 0.75]*size(SuperResImg,1))
h_axis= [0.3 0.7 0.35 0.75]*size(SuperResImg,1);
set(gca,'Xtick',[h_axis(1),h_axis(2)],'XTickLabel',{'0.3','0.7'})
set(gca,'Ytick',[h_axis(3),h_axis(4)],'YTickLabel',{'0.35','0.75'})
% set(gca,'Xtick',[h_axis(1),h_axis(2)],'XTickLabel',{'0','1'})
% set(gca,'Ytick',[h_axis(3),h_axis(4)],'YTickLabel',{'0','1'})
%title('Super-resolution Result','fontsize',8)
set(gca,'fontsize',9)
axes('position',[0.35 0.1 0.3 0.05])
colorbar('location','South')
cm = constructColormapGray;
colormap(cm)
caxis([min(c_es) max(c_es)])  
axis off
box on


axes('position',[0.7 0.22 0.3 0.65])
for ii=1:length(c_es)
    if c_es(ii)>0
          plot(xi_es(ii,1),xi_es(ii,2),'o','markersize',6,'MarkerFaceColor',[0.7 0.1 0.1])%,'MarkerEdgeColor','none')
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
axis([0.3 0.7 0.35 0.75])
%title('Peak Locations','FontSize',8)
axis square
axis xy
box on
axis on

set(gca,'fontsize',9)
%axis off
%axes('position',[0.05 0.15 0.40 0.65])
axes('position',[0.80 0.041 0.12 0.05])
plot(0.7,0,'x','color',[0.7 0.091 0.1]);
hold on
plot(0.7,-0.2,'o','markersize',6,'MarkerFaceColor',[0.7 0.1 0.1])
hold on
plot(1.5,0,'+','color',[55  150 171]./255)
plot(1.5,-0.2,'o','markersize',6,'linewidth',1)


text(0.4,0.2,'Positive','FontSize',10)
text(1.3,0.2,'Negative','FontSize',10)

text(-1.1,0.0,'Original peaks','FontSize',10)

text(-1.1,-0.2,'Estimated peaks','FontSize',10)
axis off

psize = get(gcf,'PaperSize');
wd = 16;
hg =7; 
lf = (psize(1)-wd)/2;bt = (psize(2)-hg)/2;
set(gcf,'PaperPosition',[lf bt wd hg]);
print('-dpng','-r600',['4_Sources_Pos_Neg_NoNoise_Final'])
print('-depsc2','-r700',['4_Sources_Pos_Neg_NoNoise_sameScale'])


%% Similary measure!
[SLE,Flow] = SigSimularityEMD(c_es,c(:),xi_es,[xi yi]);

%Results = [c_es xi_es];
%[num_idens,num_clusters,LE,SE,DNP,Loc_x,Loc_y,Amp] = simul_eval(Results,[xi' yi'],c',1,0.001);
