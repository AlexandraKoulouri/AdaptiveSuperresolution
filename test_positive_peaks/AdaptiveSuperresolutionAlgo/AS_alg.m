function [xi_es,c_es]= AS_alg(ob,xm,ym,sigma,alpha,grid_size,mindist,QuietPlot,xi,yi)
%Inputs
%ob(:) : measurements, low resoluton observations (inserted as 1D)
%xm(:),ym(:) coordinated of the measurements
%sigma =[sigma1 sigma2] is the standard deviations of the joint Gaussian
%convolution kernel: joint Gaussian G with G=alpha(1) G_1(sigma(1))+alpha(2) G_2(sigma(2)))
%alpha =[alpha1 alpha2] are the weights of the joint Gaussian Kernel.
%grid_size = N initial NxN uniform computational grid to solve the L1 norm problem
% MinDist = is the minimum distance between two grid points (so that the computational domain is updated) 
%QuietPlot = 1 to show the evolution of the algorithm

%Outputs:
%xi_es : estimated locations
%c_es : estimated amplitudes

%This code was created by A. Koulouri 13.5.2020

%Set the initial computational grid (uniform) in a [0 1] domain
u = linspace(0,1,grid_size);
[x,y] = meshgrid(u,u);
x0 = x(:);
y0 = y(:);
clear x y

% Convolution matrix
A = GenerateGaussianMat([xm ym],[x0 y0],sigma,alpha);

c_hat = hal(A,ob,0.01); % solve adaptive Lasso problem
%c_hat =l1_ls_nonneg(A,ob,0.01*find_lambdamax_l1_ls(A',ob))
ind_nz= find(abs(c_hat)>max(abs(c_hat))*1e-4); %thresholding - remove lower values



H = cell(1,2);
H0 = delaunay(x0,y0);

if QuietPlot ==0
    figure
    set(gcf, 'Units','centimeters', 'Position',[5 5 9 7])
    triplot(H0,x0,y0,'color',[160 160 160]./255)
    hold on
    plot(xi,yi,'xr','markersize',8)
    hold on
    plot(x0(ind_nz),y0(ind_nz),'ob','markersize',4)
    box on
    axis square
end
H{1,1} = H0;   %shows the elements
H{1,2} = 1;    %1 if we can add extra nodes, if H{i,2}=0 then we do not add extra nodes
[ind_nz_cell] = fun_arrange_nzentries_subregions(1, {length(x0)},ind_nz,c_hat);
    
%Update the support based on the non-zero entries
[x, y, updatetrue, H]=fun_update_domain({ind_nz},{x0},{y0},mindist,H); %return updatetrue=1 (if we have update the domain)
%Plot_Recursion_Result(H,x,y,{x0},{y0},ind_nz_cell,xi,yi); %plot the updated computational domain based on solution c_hat


stopiter =1;

%% Solve recursively the L1 norm problem
while updatetrue == 1      
              
      A = GenerateGaussianMat([xm ym], [cell2mat(x) cell2mat(y)],sigma,alpha);
      c_hat = hal(A,ob,0.005);
%       c_hat =l1_ls_nonneg(A,ob,0.01*find_lambdamax_l1_ls(A',ob));
%       % c_hat = lsqnonneg(A,ob);
%       yb = ob;
%       for brIter = 1:1*stopiter %increase the number of iteration as the computational grid points approach the actual points
%                 c_hat = l1_ls_nonneg(A,yb,0.02*find_lambdamax_l1_ls(A',yb));
%                 yb = yb + (ob-A*c_hat);
%                % norm(ob-A*c_hat)
%               %  if (norm(ob-A*c_hat)<= norm(ob_n-ob')); break; end
%       end
%       
      
      ind_nz = find(abs(c_hat)>max(abs(c_hat))*1e-3);
        
       
      %Create cells to store the points in the corresponding
      %sub-regions/clusters
      sizesubregion = cellfun(@size,x,'uni',false);
      nogroups = length(x);
      [ind_nz_cell] = fun_arrange_nzentries_subregions(nogroups, sizesubregion,ind_nz,c_hat);
    
      %Plot 
     % Plot_Recursion_Result(H,x,y,x,y,ind_nz_cell,xi,yi);
     
      
      [xnew, ynew, updatetrue, Hnew] = fun_update_domain(ind_nz_cell,x,y,mindist,H);  
       
     
      
      if updatetrue ==0 
          break; 
      else
           if QuietPlot==0
                Plot_Recursion_Result(Hnew,xnew,ynew,x,y,ind_nz_cell,xi,yi); %plot the updated computational domain based on solution c_hat
           end
             x = xnew; y = ynew; H = Hnew; 
      end
      
     
    stopiter =stopiter +1;
    if stopiter==30
     break;
    end
    
    
    
end


%% approximate the locations of the peaks

[ind_nz_cell,c_hat_cell]= fun_arrange_nzentries_subregions(nogroups, sizesubregion,ind_nz,c_hat); %Remove first possible empty cells
[x_f,y_f,ind_nz_cell,c_hat_cell] = remove_empty_cells(x,y,ind_nz_cell,c_hat_cell);
xi_es = zeros(length(x_f),2);
for i = 1 : length(x_f)
    xi_es(i,:) = fun_loc_est(c_hat_cell{i},x_f{i}(ind_nz_cell{i}),y_f{i}(ind_nz_cell{i}),sigma,alpha,0); %xi_es (estimated peak locations)
end

%Update the convolution matrix and estimate the amplitudes based on the
%estimated location xi_es
A = GenerateGaussianMat([xm ym], [xi_es(:,1) xi_es(:,2)],sigma,alpha);
c_es = hal(A,ob,0.01);
%c_es = lsqnonneg(A,ob);

%remove possible zero coefficients
IndZc =find(c_es==0);
xi_es(IndZc,:)=[]; 
c_es(IndZc)=[];



