function [x_est,c_est]=SDP_solver(y_f,fc,y,p,lambda)
% Script to perform super-resolution of spikes from low frequency
% measurements. More specifically, the locations of the nonzeros components 
% of x are determined from measurements of the form y=Fx where F is the low 
% pass operator that maps a function in the unit interval to its 2fc+1 lower
% Fourier series coefficients. Recovery is carried out by solving:
% 
% max_u Re(y'u) - \delta || u ||_2
% subject to max_ j |( F*u )_j | <= 1
% 
% using an equivalent SDP formulation. The script uses CVX
% (http://cvxr.com/cvx/).
%
% For more information see the paper "Super-resolution from noisy data" 
% by E. Candes and C. Fernandez-Granda. 
% Author: Carlos Fernandez-Granda 
% Email: cfgranda@stanford.edu

n=2*fc+1;
%delta_noise = sqrt(n + sqrt(2*n))*sigma; % good approximation of norm(z)
%lambda= 1*sigma*sqrt(n*log(n));
% Solve SDP
cvx_begin
%cvx_solver sdpt3
cvx_precision best
cvx_begin sdp quiet
    variable X(n+1,n+1) hermitian;
    X >= 0;
    X(n+1,n+1) == 1;
    trace(X) == 2;
    for j = 1:n-1,
        sum(diag(X,j)) == X(n+1-j,n+1);
    end
    maximize(real(X(1:n,n+1)'*y_f)-norm(X(1:n,n+1))*lambda)
cvx_end
dual_op=cvx_optval;
u = X(1:n,n+1);
aux_u =- conv(u,flipud(conj(u)));

aux_u(2*fc+1)=1+aux_u(2*fc+1);
roots_pol = roots(flipud(aux_u));
% Isolate roots on the unit circle
roots_detected = roots_pol(abs(1-abs(roots_pol)) < 1e-4);
[auxsort,ind_sort]=sort(real(roots_detected));
roots_detected = roots_detected(ind_sort);
% Roots are double so take 1 out of 2 and compute argument
t_rec_roots = angle(roots_detected(1:2:end))/2/pi;
% Argument is between -1/2 and 1/2 so convert angles
t_rec_roots(t_rec_roots < 0)= t_rec_roots(t_rec_roots < 0) + 1;  
t_rec_roots=sort(t_rec_roots);


F_est = (exp(-1i*2*pi*(-fc:1:fc)'*t_rec_roots')./sqrt(n)); % estimated Fourier matrix
length_est=length(t_rec_roots);
cvx_solver sedumi
cvx_precision best
cvx_begin quiet
    variable x_est_aux(length_est,1) complex;
    minimize(norm(x_est_aux,1))
    subject to 
        norm(y_f-sqrt(n)*F_est*x_est_aux,2)<=lambda
cvx_end
x_est_aux = sqrt(n)*x_est_aux;
primal_op=cvx_optval;
duality_gap = primal_op - dual_op;
thresh_support = 1e-4;
x_est = t_rec_roots(abs(x_est_aux)>thresh_support);
% debias estimate
%F_est_debias = exp(-1i*2*pi*x_est*(-fc:fc)); % estimated Fourier matrix
% c_est = F_est'\y;
% c_est = real(c_est);


A = DirichletKernel(fc,x_est',p);
c_est =  l1_ls(A,y,lambda);



% s = sign(real(F_est_debias*u));
% %estimate the amplitudes
% c_est = real(F_est_debias'\y - lambda*pinv(F_est_debias*F_est_debias')*s );
% % solve primal problem to check duality gap and primal support

% F_est = (exp(-1i*2*pi*(-fc:1:fc)'*x_est')); % estimated Fourier matrix


