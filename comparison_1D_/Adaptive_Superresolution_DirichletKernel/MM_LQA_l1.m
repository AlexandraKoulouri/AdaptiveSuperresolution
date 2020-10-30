function x= MM_LQA_l1(A,y,lambda,x0);

%Implementation of the L1 norm MMLGA algorithm see
%see paper for further details: 
% Comparative study of computational algorithms
%for the Lasso with high-dimensional, highly correlated data
%by Baekjin Kim1 · Donghyeon Yu2 · Joong-HoWon1 Appl Intell 2016

%created by A. Koulouri 2/4/2020


x=zeros(size(A,2),1)+x0;

N = length(x);

r=y; %residual 

At = A';

x_diff=1e20;



AtA = 2*(At*A);
Aty = 2*(At*y);

diagAtA = diag(AtA);

k=1;
maxiter = 1000;
errorTol = 1e-14;

epsilon=1e-7;
reltol = 1e-8;
pobj  = Inf; dobj  =-Inf;

quiet = 1;

%PCG parameters
pcgtol = 1e-6;
pcg_iter_max = 100;
eta = 1e-3;


while((k < maxiter) && (x_diff) > errorTol  )
    
%     for i = 1:N
%        
%         
%         
%     end
    x_old = x;
    d = 0.5*lambda./(abs(x)+epsilon);
 
    
   %M = AtA+diag(d);
    
   % P = diag(diag(M));
  
    p = diagAtA+d;
   %
     [x, pflg,prelres,pitr,presvec]= pcg(@AXfunc_l1_ls,Aty,pcgtol,pcg_iter_max,@Mfunc_l1_ls,[],x,AtA,d,1./p);
     
  %  [x, pflg,prelres,pitr,presvec] = pcg(M,Aty,pcgtol,pcg_iter_max,P,[],x_old);
    
    %do line search or something else
%     g_old = Aty-(AtA+diag(d))*x_old;
%     g  =    Aty-(AtA+diag(d))*x;
%     
%     
%     s = ((x-x_old)'*(g-g_old))/((x-x_old)'*(x-x_old));
    
    
        
    %primal and dual problem
    
     z = A*x-y;
    
    %------------------------------------------------------------
    %       CALCULATE DUALITY GAP
    %------------------------------------------------------------

    nu = 2*z;

    maxAnu = norm(At*nu./lambda,inf);
    if (maxAnu > 1)
        nu = nu/maxAnu;
    end
    pobj  =  z'*z+lambda*norm(x,1);
    dobj  =  max(-0.25*(nu'*nu)-nu'*y,dobj);
    gap   =  pobj - dobj;
    
    %Duality Gap
     if (gap/dobj < reltol) 
                    break;
     end
    
    
    
    
      if (~quiet) disp(sprintf('%4d %12.2e %15.5e %15.5e  %8d', k, gap, pobj, dobj,  pitr)); end
    
    
    x_diff=norm(x-x_old,2)/(norm(x(:),2)+eps);
    k=k+1;
    
end

%------------------------------------------------------------
%       COMPUTE AX (PCG)
%------------------------------------------------------------
function [y] = AXfunc_l1_ls(x,AtA,d,p)
%
% y = hessphi*[dx],
%
% where hessphi = [A'*A*2+D]

y = AtA*x+d.*x;

%------------------------------------------------------------
%       COMPUTE P^{-1}DX (PCG)
%------------------------------------------------------------
function [y] = Mfunc_l1_ls(x,AtA,d,p)
%
% y = P^{-1}*dx,
%

y = x.*p;

