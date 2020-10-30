%Solve the L1 norm problem is a fixed grid
close all; clear


%Peak locations and amplitudes 
xp = [4.4 13.1 10.6 5]*0.05;
yp = [2 3.2 17 8]*0.05;

c  = [1 1 1 1];

u = linspace(0,1,20);
[x,y] = meshgrid(u,u);
x = x(:);
y = y(:);


%observations
sigma = 2.5*u(2);
ob =  SimulateObservations2D(x,y,sigma,c,xp,yp);
ObImg = reshape(ob,sqrt(size(ob,2)),sqrt(size(ob,2)));
%deconvolution
% a. Basis function coefficients
A = GenerateBasisFunMat2D([x y],[x y],sigma);

% b. max regularization parameter
alpha_max = 2*max(A'*ob');

%Numerical estimation
alpha = 0.001*alpha_max;
c_hat = l1_ls_nonneg(A,ob',alpha,1e-5);

NzInd = find(c_hat>0.001*max(c_hat));


%% Results, there is not noise

figure;
set(gcf, 'Units','centimeters', 'Position',[1 1 14 7]) 
axes('position',[0.14 0.25 0.28 0.65])

%plot observations
%ObImg = reshape(ob,img_size,img_size);
%imagesc(ObImg)
colMap = colormap(hot);
colormap(colMap)
imagesc(ObImg)
title('Observations G*\mu','FontSize',8)
h_axis = axis;
set(gca,'Xtick',[h_axis(1),h_axis(2)],'XTickLabel',{'0','1'})
set(gca,'Ytick',[h_axis(3),h_axis(4)],'YTickLabel',{'0','1'})
axis square
axis on
axis xy
set(gca,'fontsize',8)
axes('position',[0.1 0.1 0.3 0.05])
colorbar('location','South')
cm = colMap;
colormap(cm)
caxis([min(ob) max(ob)])  
axis off


axes('position',[0.5 0.25 0.28 0.65])
plot(xp,yp,'xr','markersize',6.5,'linewidth',1.3);%,'color',[0.7 0.091 0.1]
hold on
u = linspace(0,1,20);
[x,y] = meshgrid(u,u);
x = x(:);
y = y(:);
plot(x,y,'.','markersize',3.7,'color',[30 89 180]/255,'markerface',[30 89 180]/255)
plot(x(NzInd),y(NzInd),'o','markersize',3,'linewidth',1);
plot(x(NzInd),y(NzInd),'.','color',[0.7 0.091 0.8])
set(gca,'fontsize',8)
title('Numerical Solution','fontsize',8)


axis square

%explanation 
axes('position',[0.55 0.001 0.32 0.23])
plot(0.7,0,'xr','linewidth',1.3)%,'color',[0.7 0.091 0.1]);
hold on
%plot(0.7,-0.2,'o','markersize',6,'MarkerFaceColor',[0.7 0.1 0.1])
plot(0.7,-0.2,'.','markersize',5.5,'color',[30 89 180]/255,'markerface',[30 89 180]/255)
hold on


plot(0.7,-0.4,'o','markersize',3,'linewidth',1)
plot(0.7,-0.4,'.','color',[0.7 0.091 0.8])
plot(0.7,-0.6,'.','color',[0.1 0.3 0.9],'linewidth',1.2,'markersize',6)

text(-1.1,0.0,'Original peaks','FontSize',8)
text(-1.1,-0.2,'Grid Points','FontSize',8)
text(-1.1,-0.4,'Non-zero coef. of \mu^N','FontSize',8)

axis off





