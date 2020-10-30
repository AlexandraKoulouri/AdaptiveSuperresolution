%In this script there are more than one peak between two grid points
%we try to recover the exact locations using over-sampling

close all
clear;

%defind interval (xmin , xmax);
xmin = 0;
xmax = 1;

%Discrete domain
Num = 35;                       % number of discrete point
Res = abs(xmax-xmin)/Num;       % resolution
p   = linspace(xmin,xmax,Num);    % values at the discrete points
x   = linspace(xmin,xmax,1000);

% std of the basis functions
sigma = 1.2*Res;

%% Generate observations
% Generate a signal with a single peak in a location xi and another one is
% the next interval
k = 15;
%locations
mu_i_1 = [p(k)+0.22* Res];
mu_i_2 = [p(k)+0.61* Res];
%amplitude
c = [1];
%observations
[ob_1,xi_1,h0] = SimulateObservations(x,sigma,mu_i_1,c);
%observations
[ob_2,xi_2,h0] = SimulateObservations(x,sigma,mu_i_2,c);


%% Generate matrix with Gaussian basis functions on the grid points and solve the munimization problem
A = GenerateBasisFunMat(x,p,sigma); 

uu = .01;
[lambda_max_1] = find_lambdamax_l1_ls_nonneg(A',ob_1);
lambda_1       = uu*lambda_max_1;

[lambda_max_2] = find_lambdamax_l1_ls_nonneg(A',ob_2);
lambda_2       = uu*lambda_max_2;

c_hat_1    = l1_ls_nonneg(A,ob_1,lambda_1,1e-6);
c_hat_2    = l1_ls_nonneg(A,ob_2,lambda_2,1e-6);



%% Plot optimality condition
%|aG(x-x1)+bG(x-x2)-G(x-x0)|<lambda
sigma2 =sqrt(2)*sigma; 
%A = GenerateBasisFunMat(linspace(xmin,xmax,50),linspace(xmin,xmax,50),sigma); 
mid_row = round(size(A,1)/2);
cc = conv(A(mid_row,:)./max(A(mid_row,:)),A(mid_row,:)./max(A(mid_row,:)));
aa = max(cc);

B =A'*A;

%just to check the convolution between two Gaussian distributions
[sigma1,mu1,aa]=mygaussfit(p,B(:,round(size(B,2)*0.5)));
G = @(u)aa* exp(-0.5 * ((u)./sigma2).^2);%./ (sqrt(2*pi) .* sigma2);
x_int = linspace(p(k-2),p(k+3),1400); %define a small interval between the actual peak

%%

F_1 = 0; % Optimality condition values around the actual peak location
F_2 = 0;
for j = 1:length(xi_1)
   F_1 = c(j)*G(x_int-xi_1(j));
end

for j = 1:length(c_hat_1);
   F_1 = F_1 -c_hat_1(j) *G(x_int-p(j)); 
end

for j = 1:length(xi_2)
   F_2 = c(j)*G(x_int-xi_2(j));
end

for j = 1:length(c_hat_1);
   F_2 = F_2 -c_hat_2(j) *G(x_int-p(j)); 
end

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

axes('position',[0.07 0.15 0.38 0.57])
plot(x_int,F_1*2/lambda_1-1,'color',[0.1 0.5 0.7],'linewidth',0.6); hold on
stem(p(k-2:k+2),zeros(length(p(k-2:k+2)),1),'o','markersize',2.5,'color',[30 89 180]/255,'markerface',[30 89 180]/255);
hold on

plot(xi_1,0,'x','color',[0.72 0.12 0.1],'markersize',8.5)
plot(p(k),0,'o','markersize',6.5,'color',[0.1 0.3 0.9])
%title('\mu^N with one non-zero coefficients','FontSize',7)
axes('position',[0.25 0.07 0.001 0.001])
text(-1.5,0.425,'(A)')

axes('position',[0.55 0.15 0.38 0.57])
plot(x_int,F_2*2/lambda_2-1,'color',[0.1 0.5 0.7],'linewidth',0.6); hold on
stem(p(k-2:k+2),zeros(length(p(k-2:k+2)),1),'o','markersize',2.5,'color',[30 89 180]/255,'markerface',[30 89 180]/255);
hold on

plot([p(k), p(k+1)],[0 0],'o','markersize',6.5,'color',[0.1 0.3 0.9])
plot(xi_2,0,'x','color',[0.72 0.12 0.1],'markersize',8.5)
%title('\mu^N with two non-zero coefficients','FontSize',7)
axes('position',[0.73 0.07 0.001 0.001])
text(-1.5,0.825,'(B)')





% set(gcf,'PaperUnits','centimeters')
% psize = get(gcf,'PaperSize');
% wd = 14;hg =8; lf = (psize(1)-wd)/2;bt = (psize(2)-hg)/2;
% set(gcf,'PaperPosition',[lf bt wd hg]);
% print('-depsc2','-r300',['OptimalityCurves'])