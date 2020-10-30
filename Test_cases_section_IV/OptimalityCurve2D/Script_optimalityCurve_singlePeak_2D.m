%Script: plot the optimality curve
clear
addpath(genpath(cd)) 

close all; clear

%Image
img_size = 30;

u = linspace(0,1,img_size);
hsize = u(2)-u(1);

[x,y] = meshgrid(u,u);
x = x(:);
y = y(:);

%remove some nodes
%x([435 436 465 466])=[];
%y([435 436 465 466])=[];
% x(NzInd)=[];
% y(NZInd)=[];
%Peak locations and amplitudes 
xp = [y(round(img_size/2))+0.2*hsize]';
yp = [y(round(img_size/2))+0.22*hsize]';

c  = [1];

%observations
sigma = 1.7*hsize;
ob =  SimulateObservations2D(x,y,sigma,c,xp,yp);

%ObImg = reshape(ob,sqrt(size(ob,2)),sqrt(size(ob,2)));
%deconvolution
% a. Basis function coefficients
A = GenerateBasisFunMat2D([x y],[x y],sigma);



%% Numerical estimation
uu = .01;
[lambda_max] = 2*max(A'*ob');
lambda       = uu*lambda_max;
residual = 1;
yy =ob';
c_hat    = l1_ls_nonneg(A,yy,lambda,1e-6);
NzInd = find(c_hat>1e-4*max(c_hat));

%% Optimality curve 2D F(x) = sum c_i H(x-xi) - sum c_hat H(x-x_i)

C = reshape(A(:,651),img_size,img_size);
aa = max(max(C'*C));  
Sigma = [sigma 0;0 sigma];
Sigma1 = Sigma;
Sigma2 = Sigma;
Sigma_c =[2*sigma^2 0;0 2*sigma^2]^(-1);
H = @(u)aa^2*exp(-0.5*u'*Sigma_c*u); 
res = 75;
x_int = linspace(xp-2.5*hsize,xp+3*hsize,res);
y_int = linspace(yp-2.5*hsize,yp+3*hsize,res);
[xint1,yint1] = meshgrid(x_int,y_int);
xint = xint1(:);
yint = yint1(:);

F = zeros(size(xint,1),1);
U = 0;
for i = 1: size(xint,1)
    for j = 1:length(c)
        U = U+ c(j)*H([xint(i)-xp(j) yint(i)-yp(j)]');
    end
    F(i) = U;
      U = 0;
%F(i) = c*H([xint(i)-xp yint(i)-yp]');
end 
U = 0;
%F2=zeros(size(xint,1),1);
for i = 1:size(xint,1);
    for jj = 1:length(c_hat)
    
        U = U + c_hat(jj)*H([xint(i)-x(jj) yint(i)-y(jj)]');
    end
   F(i) = F(i) - U;
   U = 0;
end

F = 2*F/lambda-1;


  %% 
close all



%%
figure;%('position',[25 290 900 800]);
set(gcf, 'Units','centimeters', 'Position',[5 5 14 8]) 

axes('position',[0.01 0.7 0.1 0.25])
plot(0.4,0.0,'o','markersize',1.5,'color',[30 89 180]/255,'markerface',[30 89 180]/255)
hold on
text(1.4,-0.2,'Non-zero coef. of \mu^N','FontSize',8)
plot([-.1 1.],[-.4 -.4], 'color',[0.1 0.5 0.7]);hold on
plot(0.4,-0.6,'x','markersize',8.5,'color',[0.72 0.12 0.1])
plot(0.4, -0.2,'o','markersize',5.5,'color',[0.1 0.3 0.9])
text(1.4,0,['Grid points'],'FontSize',8)
text(1.5,-.4,'p(x), (p"(x)<0)','FontSize',8)
text(1.4,-0.6,'\xi','FontSize',8)

axis off
%
% 
axes('position',[0.07 0.15 0.38 0.57])



%load colormapInvGray
gray = constructColormapGray;
FImg = reshape(F,res,res);
h=surf(yint1,xint1,FImg,'EdgeColor','none')
colormap gray
%beta = .3;
alpha(0.1)
%brighten(beta);


hold on
xs = linspace(xp-2.5*hsize,xp+3*hsize,10);
ys = linspace(yp-2.5*hsize,yp+3*hsize,10);
[x_s y_s]= meshgrid(xs,ys);

surf(x_s,y_s,zeros(size(x_s))-0.02,'EdgeColor','none')
%alpha(0.3)
plot3(yp,xp,max(max(FImg))+0.1,'x','color',[0.72 0.12 0.1],'markersize',10.5)

plot3(x(NzInd),y(NzInd),zeros(length(NzInd),1)+0.1,'o','markersize',6.5,'color',[0.1 0.3 0.9])%plot non-zero points
dist1= (x-x(NzInd(1))).^2+(y-y(NzInd(1))).^2;
indkeep = find(dist1<=0.02*max(dist1));


plot3(x(indkeep),y(indkeep),zeros(length(x(indkeep)),1)+0.1,'.','markersize',6.5,'color',[30 89 180]/255,'markerface',[30 89 180]/255);

% hold on
% 
% plot([p(k), p(k+1)],[0 0],'o','markersize',6.5,'color',[0.1 0.3 0.9])



% xlabel('d_x^{(1)}','FontSize',12)
% ylabel('d_y^{(1)}','FontSize',12)
zlabel('p(\bf{x})','FontSize',12)
%text(X(1,1),Y(1,1),0.01*F(end)+0.01,'a=0.5a_{max}','FontSize',12)
view(-33,27)
axis square
axis tight

set(get(gca,'ZLabel'),'Rotation',0);
set(gcf, 'renderer', 'painters');
set(gca,'ZTick',[0])

axes('position',[0.55 0.15 0.38 0.57])
 
[C,h]=contour(x_int,y_int,FImg);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*3)
clabel(C,h);
hold on
plot3(yp,xp,max(max(FImg))+0.1,'x','color',[0.72 0.12 0.1],'markersize',10.5)

plot3(x(NzInd),y(NzInd),zeros(length(NzInd),1)+0.1,'o','markersize',6.5,'color',[0.1 0.3 0.9])%plot non-zero points
dist1= (x-x(NzInd(1))).^2+(y-y(NzInd(1))).^2;
indkeep = find(dist1<=0.02*max(dist1));


plot3(x(indkeep),y(indkeep),zeros(length(x(indkeep)),1)+0.1,'.','markersize',6.5,'color',[30 89 180]/255,'markerface',[30 89 180]/255);


set(gcf,'PaperUnits','centimeters')
psize = get(gcf,'PaperSize');
wd = 14;hg =8; lf = (psize(1)-wd)/2;bt = (psize(2)-hg)/2;
set(gcf,'PaperPosition',[lf bt wd hg]);
print('-depsc2','-r300',['OptimalityCurves2D'])

%save figures
% print(gcf, '-dpdf', 'ShapeCostFun.pdf');
% print(gcf, '-dpng', 'ShapeCostFun.png');
% print(gcf, '-depsc2', 'ShapeCostFun.eps');



