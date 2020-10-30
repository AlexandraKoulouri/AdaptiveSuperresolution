function [xi_es,c_es] = AS_1D(ob,grid_size,z,fc,xi,c,Lim,alpha)


%ob: observations
%grid_size: original computational grid
%p: measurement points
%true location xi and amplitude c (just for visualization)
%code for Dirichlet kernel

%In this code we iteratively update the grid domain around the locations
%were there non-zero entries in solution c_hat. In 1D implementation we do
%not use elements and nodes as in the 2D implementation(but that would be useful to
%optimize the clustering process)

%return: locations xi_es and amplitudes c_es

% first computational grid
x = linspace(min(z)+0.0001,max(z)-0.001,grid_size);


%% Generate matrix with Gaussian basis functions on the grid points and solve the munimization problem
A = DirichletKernel(fc,x,z);

 [lambda_max] = find_lambdamax_l1_ls(A',ob);
 lambda       = alpha*lambda_max;

c_hat    = l1_ls(A,ob,lambda,1e-6);
%c_hat = hal(A,ob,alpha);

Ind_nz = find(abs(c_hat)>=1e-4*max(c_hat));
 
% figure;
% set(gcf, 'Units','centimeters', 'Position',[5 5 12 6])
% stem(xi,c,'rx','linewidth',1.2,'markersize',6.5,'linestyle','-.'); hold on
% %stem(z,zeros(length(z,1),1),'.','markersize',9.5)
% plot(z,ob,'color',[0.1 0.5 0.7],'linewidth',0.6)
% hold on
% stem(x(Ind_nz),c_hat(Ind_nz),'color',[0.1 0.5 0.8],'markersize',5,'linewidth',1.2);



%% Add extra points
%x is a vector whth the total number of points (the actual grid points and the added ones)
%x_ls includes only the candidate locations where the coefficients of the
%de-convolution problem will be estimated in the next step
[x, ~ ,x_ls] = AddExtraPoints_1D(Ind_nz,x,Lim);

x_ls_cell = cell(3,1);
c_hat_ls_cell = cell(3,1);
Ind_nz_cell = cell(3,1);

bool = 0;
iter = 0;
 while bool~=1
         iter = iter+1;
        
         %estimate basis functions
         A = DirichletKernel(fc,x_ls',z);
         %estimate least squares solution 
         [lambda_max] = find_lambdamax_l1_ls(A',ob);
         lambda   = alpha*lambda_max;
        c_hat   = l1_ls(A,ob,lambda,1e-6);   
       % c_hat = hal(A,ob,0.1*alpha);

        %c_hat = l1_ls_nonneg(A,ob,lambda,1e-6);   
        
         %c_hat_LS = lsqlin(A_LS,ob_n,-eye(length(x_ls)),zeros(length(x_ls),1),[],[],[],[],[],'LargeScale','off');
         Ind_nz = find(abs(c_hat)>0.0001*max(c_hat)); % indices of the non zero coefficients
         
         x_ls_cell(iter) = {x_ls};
         c_hat_ls_cell(iter)  = {c_hat};
         Ind_nz_cell(iter) = {Ind_nz};
         
% 
%           figure;
%           set(gcf, 'Units','centimeters', 'Position',[5 5 12 6])
%           stem(xi,c,'rx','linewidth',1,'markersize',5.5,'linestyle','-.')
%           hold on
%           plot(z,ob,'color',[0.1 0.5 0.7],'linewidth',0.6)
%           stem(x_ls(Ind_nz), c_hat(Ind_nz),'o','markersize',5,'color',[0 0 0])
%          
                       
          x_ls_t = x_ls(Ind_nz);
          Ind_ls = find(ismember(x,x_ls_t) ==1);
          [x_new,bool,x_ls_t] = AddExtraPoints_1D(Ind_ls,x,Lim);
        % plot(x_ls_t,zeros(length(x_ls_t),1),'x')
          if bool == 0; x = x_new; x_ls = x_ls_t; end
 end
 
   
 %% Estimate the peak locations
 %estimate the amplitude of the convoluted basis functions
% Ind_nz = find(abs(c_hat)>0.0001*max(c_hat)); % indices of the non zero coefficients
%   x_ls_t = x_ls(Ind_nz);       
 
 mid_row = round(size(A,1)/2);
 cc = conv(A(mid_row,:),A(mid_row,:));
 [~,xi_es] = fun_EstPeakLocation(x_ls,c_hat,max(cc),lambda);
 
 A = DirichletKernel(fc,xi_es,z);
 %c_es = lsqlin(A,ob,-eye(length(xi_es)),zeros(length(xi_es),1),[],[],[],[],[],'LargeScale','off');
%[lambda_max] = find_lambdamax_l1_ls(A',ob);
% lambda   = 0.0001*lambda_max;
 %c_es   = l1_ls(A,ob,lambda,1e-6);    %remove noise bias run bregman iteration - option
%c_es= lsqnonneg(A,ob); %when positivity is applied
%c_es = hal(A,ob,alpha);
c_es = A\ob;
ind = find(abs(c_es)<1e-3);
c_es(ind)=[];
xi_es(ind)=[];

% figure;
% stem(xi,c,'rx','linewidth',1,'markersize',5.5,'linestyle','-.')
% hold on
% plot(z,ob,'color',[0.1 0.5 0.7],'linewidth',0.6)
%   
% hold on
% stem(xi_es, c_es,'ok','markersize',5,'color',[0 0 0])

% [xi,ii]= sort(xi);
% xi_es1 = sort(xi_es);
 
 
 