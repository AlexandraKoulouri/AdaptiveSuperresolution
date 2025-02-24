function c_hat = hal(A,b,lambda);

%A = GenerateBasisFunMat2D([xm ym], [cell2mat(x) cell2mat(y)],sigma);

%this code is based on the Hierarchical adaptive lasso (see Murphy's book:
%Machine Learning: A Probabilistic Perspective)

%this code was created by A. Koulouri 5.5.2019, updated in 24.2.2025

%lambda takes values in (0,1) interval

lambda =lambda*find_lambdamax_l1_ls(A',b); %adjust the regularization parameter (maybe in the future!)
beta = 0.0001;
kappa = lambda*beta-1;

c_hat_old = zeros(size(A,2),1)+1e6;
stop_cond = false;
cost_old = norm(A*c_hat_old-b,2)^2+kappa*sum(log(abs(c_hat_old)/beta+1));%lambda*norm(c_hat_old,1);

iter = 1; max_iter=20;


while(~stop_cond) 
    
              
       % c_hat=l1_ls(A,b,lambda);
        c_hat= MM_LQA_l1(A,b,lambda,zeros(size(A,2),1));
        cost_new = norm(A*c_hat-b,2)^2+kappa*sum(log(abs(c_hat)/beta+1)); %norm(lambda.*c_hat,1);
        lambda = (kappa+1)./(beta+abs(c_hat));
       
        delta_cost = abs(cost_old-cost_new);
        avg_cost = (abs(cost_new)+abs(cost_old)+eps)/2;
        if (delta_cost/avg_cost)< 1e-4
        stop_cond= true;
        end
%         if  (cost_old-cost_new)<-2*eps
%         warning('convergenceTest:neglog_Decrease', 'objective increased!'); 
%         end    
        cost_old = cost_new;
        c_hat_old=c_hat;
        iter = iter+1;
        if iter>max_iter
            break;
        end
end
