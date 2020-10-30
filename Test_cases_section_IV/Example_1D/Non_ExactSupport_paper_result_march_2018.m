%In this script there are more than one peak between two grid points
%we try to recover the exact locations 
close all
clear;

%defind interval (xmin , xmax);
xmin = 0;
xmax = 1;

%Discrete domain
Num = 51;                       % size of computational grid
Res = abs(xmax-xmin)/Num;       % resolution
x = linspace(xmin,xmax,Num+1);  % grid points

p   = linspace(xmin,xmax,100);  % observationpoints

% std of the basis functions
sigma = 3*p(2);

%True signal ->
%locations
mu_i = [0.23 0.58 0.83];
%amplitude
c = [0.5 0.9 0.7 ];
%Simulated observations
[ob,xi,h0] = SimulateObservations(p,sigma,mu_i,c);

%% Generate matrix with Gaussian basis functions on the grid points and solve the munimization problem
A = GenerateBasisFunMat(p,x,sigma); 

[lambda_max] = find_lambdamax_l1_ls_nonneg(A',ob);
lambda       = 0.01*lambda_max;

residual = 1;
y =ob;
c_hat    = l1_ls_nonneg(A,y,lambda,1e-6); %use positivity here
Ind = find(c_hat>=0.0015*max(c_hat));


%% Plot numerical solution

figure;%('position',[25 290 900 800]);
set(gcf, 'Units','centimeters', 'Position',[5 5 14 7])

axes('position',[0.045 0.2 0.27 0.55])
%texthvc(0,0,['Resolution = ', num2str(Res)], 'cl',[0 0 0])
stem(xi,c,'x','linewidth',.65,'markersize',4.5,'color',[0.72 0.12 0.1])%,'linestyle','-.');
hold on
stem(x,zeros(length(x),1),'o','markersize',2.5,'color',[30 89 180]/255,'markerface',[30 89 180]/255)%[0,110,0]./255)
title('Test Case: Non-Exact Support','FontSize',8)
%axis([p(k-2) p(k+2) 0 1.4] )

axes('position',[0.05 0.70 0.015 0.4])
plot(0,0.0,'.','markersize',8.5,'color',[30 89 180]/255)
%texthvc(0,0,['Grid points (N=',num2str(Num),')'], 'cl',[0 0 0])
text(1.2,0,['Grid points (N=',num2str(Num+1),')'],'FontSize',8)
axis off



 axes('position',[0.376 0.2 0.27 0.55])
 plot(p,ob,'color',[0.1 0.5 0.7],'linewidth',.65)%,'color',[0,110,0]./255)
 hold on
% stem(p,zeros(length(p)),'.','markersize',7.5)
 %stem(xi,c,'rx','linewidth',1.1,'markersize',4.5,'linestyle','-.')
 title('Observations G*\mu','FontSize',8)
 %set(gca,'Ytick',0,'YTickLabel','')
axis(gca)

 
 axes('position',[0.72 0.2 0.27 0.55])
 stem(x(Ind),c_hat(Ind),'color', [0.1 0.3 0.9],'markersize',4.5,'linewidth',.65);
 hold on
 %stem(xi,c,'rx','linewidth',1,'markersize',7.5,'linestyle','-.')
stem(x,zeros(length(x),1),'o','markersize',2.5,'color',[30 89 180]/255,'markerface',[30 89 180]/255)
 %stem(p(k-2:k+2),zeros(length(p(k-2:k+2)),1),'.','markersize',9.5);
 title(['Reconstruction \mu^N '],'FontSize',8)
 
 %% Estimate  peak locations
 %estimate the amplitude of the convoluted basis functions
mid_row = round(size(A,1)/2);
cc = conv(A(mid_row,:),A(mid_row,:));
[c_es,xi_es] = fun_EstPeakLocation(x,c_hat,sigma,lambda,max(cc));

%plot true locations and amplitudes and reconstructed ones
figure;%('position',[25 290 900 800]);
set(gcf, 'Units','centimeters', 'Position',[5 5 12 7])
 axes('position',[0.1 0.2 0.75 0.6])
% stem(p(k-2:k+2),zeros(length(p(k-2:k+2)),1),'.','markersize',9.5);
 stem(xi,c,'x','color','r','linewidth',0.8,'markersize',4.5,'linestyle','-.')
 hold on;
 stem(xi_es,c_es,'color','k','linewidth',0.8,'markersize',5); 
 %set(gca,'ylim',x_int)
 stem(x,zeros(length(x),1),'.','markersize',7.5)

 title('Peak Estimation','FontSize',8)
 [xi,ii]= sort(xi);
 xi_es1 = sort(xi_es);
 h = legend('Actual Peaks','Estimated Peaks')%'FontSize',6)
 
 set(h,'FontSize',6,'Location','northwestoutside','box','off')
 set(gcf,'PaperUnits','centimeters')
 psize = get(gcf,'PaperSize');
 wd =12;hg =7; lf = (psize(1)-wd)/2;bt = (psize(2)-hg)/2;
 set(gcf,'PaperPosition',[lf bt wd hg]);
 %print('-depsc2','-r300',['NonExactSupportCases_PeakExtraction_mix'])

